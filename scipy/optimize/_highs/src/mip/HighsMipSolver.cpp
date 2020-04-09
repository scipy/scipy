/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolver.h"

#include "lp_data/HighsModelUtils.h"

// Branch-and-bound code below here:
// Solve a mixed integer problem using branch and bound.
HighsMipStatus HighsMipSolver::runMipSolver() {
  std::cout << "Warning: HiGHS MIP solver is under construction at the moment."
            << std::endl
            << "Running HiGHS MIP solver..." << std::endl;

  // Start timer.
  timer_.startRunHighsClock();
  double mip_solve_initial_time = timer_.readRunHighsClock();

  // Load root node lp in highs and turn printing off.
  passModel(mip_);
  // Set the options for this Highs instance according to the options for the
  // MIP solver
  options_ = options_mip_;

  // Report deviations from default options settings
  writeHighsOptions("");
  options_.message_level = 0;

  const bool only_write_as_mps = false;  // true;//
  if (only_write_as_mps) {
    printf("Only writing out the MIP as MPS\n");
    writeModel("mip.mps");
    //    options_.message_level=7; printf("Writing out the MIP on stdout\n");
    //    writeModel(""); options_.message_level=4;
    return HighsMipStatus::kUnderDevelopment;
  }
  HighsMipStatus root_solve_status = solveRootNode();
  num_nodes_solved++;
  root_objective_ = info_.objective_function_value;
  reportMipSolverProgressLine("", true);
  reportMipSolverProgress(root_solve_status);
  if (root_solve_status != HighsMipStatus::kRootNodeOptimal)
    return root_solve_status;

  // Start tree by making root node.
  // Highs ignores integrality constraints.
  Node root(-1, 0.0, 0, 0);
  root.col_lower_bound = lp_.colLower_;
  root.col_upper_bound = lp_.colUpper_;
  root.integer_variables = lp_.integrality_;
  root.primal_solution = solution_.col_value;
  root.objective_value = info_.objective_function_value;

  //  writeSolutionForIntegerVariables(root);

  // Add and solve children.
  HighsMipStatus tree_solve_status = solveTree(root);
  reportMipSolverProgress(tree_solve_status);

  // Stop and read the HiGHS clock, then work out time for this call
  double mip_solve_final_time = timer_.readRunHighsClock();

  int num_nodes_formed = tree_.getNumNodesFormed();
  int num_integer_solutions = tree_.getNumIntegerSolutions();
  int num_nodes_unsolved =
      num_nodes_formed - num_nodes_solved - num_nodes_pruned;
  HighsPrintMessage(options_mip_.output, options_mip_.message_level, ML_MINIMAL,
                    "\nMIP solver summary\n");
  HighsPrintMessage(options_mip_.output, options_mip_.message_level, ML_MINIMAL,
                    "Number of nodes formed   = %9d\n", num_nodes_formed);
  HighsPrintMessage(options_mip_.output, options_mip_.message_level, ML_MINIMAL,
                    "Number of nodes solved   = %9d\n", num_nodes_solved);
  HighsPrintMessage(options_mip_.output, options_mip_.message_level, ML_MINIMAL,
                    "Number of nodes pruned   = %9d\n", num_nodes_pruned);
  HighsPrintMessage(options_mip_.output, options_mip_.message_level, ML_MINIMAL,
                    "Number of nodes unsolved = %9d\n", num_nodes_unsolved);

  HighsPrintMessage(options_mip_.output, options_mip_.message_level, ML_MINIMAL,
                    "Number of IFS found      = %9d\n", num_integer_solutions);

  if (tree_.getBestSolution().size() > 0) {
    if (num_nodes_unsolved)
      HighsPrintMessage(options_mip_.output, options_mip_.message_level,
                        ML_MINIMAL, "ERROR: number of nodes unsolved = %9d\n",
                        num_nodes_unsolved);
    hmos_[0].unscaled_model_status_ = HighsModelStatus::OPTIMAL;
    std::stringstream message;
    message << std::endl;
    message << "Run status : " << highsMipStatusToString(tree_solve_status)
            << std::endl;
    message << "Objective  : " << std::scientific << tree_.getBestObjective()
            << std::endl;
    message << "Time       : " << std::fixed << std::setprecision(3)
            << mip_solve_final_time - mip_solve_initial_time << std::endl;
    message << std::endl;

    HighsPrintMessage(options_mip_.output, options_mip_.message_level,
                      ML_MINIMAL, message.str().c_str());
    // todo: handle feasible vs optimal case once you have a timeout.
  } else {
    HighsPrintMessage(options_mip_.output, options_mip_.message_level,
                      ML_MINIMAL, "No feasible solution found.");
    // todo: handle infeasible vs timeout case once you have a timeout.
  }

  if (tree_solve_status != HighsMipStatus::kTreeExhausted) {
    std::cout << "Warning: tree not explored entirely." << std::endl;
    return tree_solve_status;
  }

  return HighsMipStatus::kUnderDevelopment;
}

#ifdef HiGHSDEV
void HighsMipSolver::writeSolutionForIntegerVariables(Node& node) {
  for (int iCol = 0; iCol < lp_.numCol_; iCol++) {
    if (!lp_.integrality_[iCol]) continue;
    printf("%2d [%10.4g, %10.4g, %10.4g]\n", iCol, node.col_lower_bound[iCol],
           node.primal_solution[iCol], node.col_upper_bound[iCol]);
  }
}
#endif

HighsMipStatus HighsMipSolver::solveNode(Node& node, bool hotstart) {
  // Force calls within run() to be silent by setting the HiGHS
  // logfile to NULL and the HiGHS message_level to zero.
  bool no_highs_log = true;
  // Use save_message_level and save_logfile to keep the original
  // message_level and logfile in case they are changed
  std::string save_presolve;
  int save_message_level;
  FILE* save_logfile;

  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsModelStatus use_model_status = HighsModelStatus::NOTSET;

  // When full_highs_log is true, run() is verbose - for debugging
  bool full_highs_log = false;
  // Setting check_node_id forces full logging for a particular node
  const int check_node_id = HIGHS_CONST_I_INF;  // 33663;//

  //     printf("SolveNode: Id = %d; ParentId = %d; BranchCol = %d\n", node.id,
  //     node.parent_id, node.branch_col);
  if (node.id == check_node_id) {
    // Switch on full logging for this node - {} used so VScode can
    // stop on this line
    full_highs_log = true;
    //    options_.mip_report_level = 2;
    printf("node%d: %d; %d\n", check_node_id, lp_.numCol_, lp_.numRow_);
  }
  if (hotstart) {
    // Apply changes to LP from node. For the moment only column bounds.
    // Get the original message_level and logfile in case they are set to
    // something different for run()
    save_message_level = options_.message_level;
    save_logfile = options_.logfile;
    save_presolve = options_.presolve;
    if (full_highs_log) {
      // Using full logging, so prevent "no logging"
      no_highs_log = false;
      options_.message_level = 7;
      options_.logfile = stdout;
    }
    if (no_highs_log) {
      // Using no logging, so prevent it
      options_.message_level = 0;
      options_.logfile = NULL;
    }

    changeColsBounds(0, mip_.numCol_ - 1, &node.col_lower_bound[0],
                     &node.col_upper_bound[0]);

    if (node.id == check_node_id) {
      printf("Writing node%1d.mps\n", check_node_id);
      writeModel("node33663.mps");
      //      basis_.valid_ = false; options_.presolve = on_string;
    }

    call_status = run();
    return_status = interpretCallStatus(call_status, return_status, "run()");
    if (return_status == HighsStatus::Error) return HighsMipStatus::kNodeError;

    call_status = getUseModelStatus(use_model_status,
                                    unscaled_primal_feasibility_tolerance,
                                    unscaled_dual_feasibility_tolerance, true);
    return_status =
        interpretCallStatus(call_status, return_status, "getUseModelStatus()");
    if (return_status == HighsStatus::Error) return HighsMipStatus::kNodeError;

    // Reset the values of message_level and logfile
    options_.message_level = save_message_level;
    options_.logfile = save_logfile;
    options_.presolve = save_presolve;
    const bool check_hotstart = false;  // true;//
    if (check_hotstart) {
      HighsModelStatus hotstart_model_status = use_model_status;
      double hotstart_objective = info_.objective_function_value;

      Highs highs;
      if (no_highs_log) {
        highs.options_.logfile = NULL;
        highs.options_.message_level = 0;
      } else if (full_highs_log) {
        highs.options_.message_level = 4;
      }
      HighsLp lp_node = mip_;
      lp_node.colLower_ = node.col_lower_bound;
      lp_node.colUpper_ = node.col_upper_bound;
      highs.passModel(lp_node);

      highs.options_.presolve = off_string;
      if (node.id == check_node_id) highs.options_.presolve = on_string;
      call_status = highs.run();
      return_status = interpretCallStatus(call_status, return_status, "run()");
      if (return_status == HighsStatus::Error)
        return HighsMipStatus::kNodeError;

      HighsModelStatus check_model_status;
      call_status = highs.getUseModelStatus(
          check_model_status, unscaled_primal_feasibility_tolerance,
          unscaled_dual_feasibility_tolerance);
      return_status = interpretCallStatus(
          call_status, return_status, "getUseModelStatus(use_model_status)");
      if (return_status == HighsStatus::Error)
        return HighsMipStatus::kNodeError;

      double check_objective = highs.info_.objective_function_value;
      if (check_model_status != hotstart_model_status) {
        // Check whether the model status is the same
        printf(
            "SolveNode ERROR: %s = check_model_status != hotstart_model_status "
            "= %s\n",
            highsModelStatusToString(check_model_status).c_str(),
            highsModelStatusToString(hotstart_model_status).c_str());
      } else if (check_model_status == HighsModelStatus::OPTIMAL) {
        // Check that the same optimal objective value is found
        double dl_objective = fabs(check_objective - hotstart_objective) /
                              max(1.0, fabs(check_objective));
        double tl_dl_objective = 1e-9;
        if (dl_objective > tl_dl_objective || full_highs_log) {
          printf("SolveNode");
          if (dl_objective > tl_dl_objective) {
            printf(" ERROR");
          }
          printf(
              ": Optimal objective difference = %g from (hotstart; check) = "
              "(%g; %g)\n",
              dl_objective, check_objective, hotstart_objective);
        }
      }
    }
  } else {
    // solve from scratch to test
    Highs highs;
    highs.options_.message_level = 0;
    HighsLp lp_node = mip_;
    lp_node.colLower_ = node.col_lower_bound;
    lp_node.colUpper_ = node.col_upper_bound;
    highs.passModel(lp_node);
    call_status = highs.run();
    return_status = interpretCallStatus(call_status, return_status, "run()");
    if (return_status == HighsStatus::Error) return HighsMipStatus::kNodeError;

    call_status = highs.getUseModelStatus(use_model_status,
                                          unscaled_primal_feasibility_tolerance,
                                          unscaled_dual_feasibility_tolerance);
    return_status = interpretCallStatus(call_status, return_status,
                                        "getUseModelStatus(use_model_status)");
    if (return_status == HighsStatus::Error) return HighsMipStatus::kNodeError;
  }

  switch (return_status) {
    case HighsStatus::Warning:
      if (use_model_status == HighsModelStatus::REACHED_TIME_LIMIT)
        return HighsMipStatus::kTimeout;
      if (use_model_status == HighsModelStatus::REACHED_ITERATION_LIMIT)
        return HighsMipStatus::kReachedSimplexIterationLimit;
      return HighsMipStatus::kNodeNotOptimal;
    case HighsStatus::Error:
      return HighsMipStatus::kNodeError;
    default:
      break;
  }

  switch (use_model_status) {
    case HighsModelStatus::OPTIMAL:
      node.primal_solution = solution_.col_value;
      node.objective_value = info_.objective_function_value;
      return HighsMipStatus::kNodeOptimal;
    case HighsModelStatus::PRIMAL_INFEASIBLE:
      return HighsMipStatus::kNodeInfeasible;
    case HighsModelStatus::PRIMAL_UNBOUNDED:
      return HighsMipStatus::kNodeUnbounded;
    case HighsModelStatus::REACHED_TIME_LIMIT:
      return HighsMipStatus::kTimeout;
    case HighsModelStatus::REACHED_ITERATION_LIMIT:
      return HighsMipStatus::kReachedSimplexIterationLimit;
    case HighsModelStatus::NOTSET:
      return HighsMipStatus::kNodeError;
    default:
      printf("HighsModelStatus: %s\n",
             highsModelStatusToString(use_model_status).c_str());
      break;
  }

  return HighsMipStatus::kNodeNotOptimal;
}

HighsMipStatus HighsMipSolver::solveRootNode() {
  HighsStatus lp_solve_status = HighsStatus::Error;
  HighsModelStatus use_model_status = HighsModelStatus::NOTSET;
  //  options_.presolve = off_string;
  bool no_highs_log = true;  // false;//
  int save_message_level = options_.message_level;
  FILE* save_logfile = options_.logfile;
  if (no_highs_log) {
    options_.logfile = NULL;
    options_.message_level = 0;
  }
  lp_solve_status = run();
  if (no_highs_log) {
    options_.logfile = save_logfile;
    options_.message_level = save_message_level;
  }
  options_.presolve = off_string;
  use_model_status = model_status_;

  switch (lp_solve_status) {
    case HighsStatus::Warning:
      return HighsMipStatus::kRootNodeNotOptimal;
    case HighsStatus::Error:
      return HighsMipStatus::kRootNodeError;
    default:
      break;
  }

  if (use_model_status != HighsModelStatus::OPTIMAL)
    return HighsMipStatus::kRootNodeNotOptimal;

  return HighsMipStatus::kRootNodeOptimal;
}

HighsMipStatus HighsMipSolver::solveTree(Node& root) {
  // The method branch(...) below calls chooseBranchingVariable(..) which
  // currently returns the first violated one. If a branching variable is found
  // children are added to the stack. If there are no more violated integrality
  // constraints we have a feasible solution, if it is best than current best,
  // the current best is updated.

  tree_.setMipReportLevel(options_.mip_report_level);

  tree_.branch(root);

  // While stack not empty.
  //   Solve node.
  //   Branch.
  while (!tree_.empty()) {
    if (timer_.readRunHighsClock() > options_.time_limit)
      return HighsMipStatus::kTimeout;
    if (tree_.getNumNodesFormed() > options_.mip_max_nodes)
      return HighsMipStatus::kMaxNodeReached;
    Node& node = tree_.next();
    double best_objective;
    best_objective = tree_.getBestObjective();
    if (node.parent_objective >= best_objective) {
      // Don't solve if we can't better the best IFS
      if (options_.mip_report_level > 1)
        printf("Don't solve since no better than best IFS of %10.4g\n",
               best_objective);
      num_nodes_pruned++;
      tree_.pop();
      continue;
    }
    HighsMipStatus node_solve_status = solveNode(node);
    num_nodes_solved++;

    switch (node_solve_status) {
      case HighsMipStatus::kNodeOptimal:
        reportMipSolverProgress(node_solve_status);
        if (options_.mip_report_level > 1) {
          printf("Node %9d (branch on %2d) optimal objective %10.4g: ", node.id,
                 node.branch_col, node.objective_value);
          /*
            std::cout << "Node " << node.id
            << " solved to optimality." << std::endl;
          */
        }
        tree_.pop();
        // Don't branch if we can't better the best IFS
        double best_objective;
        best_objective = tree_.getBestObjective();
        if (node.objective_value >= best_objective) {
          if (options_.mip_report_level > 1)
            printf("Don't branch since no better than best IFS of %10.4g\n",
                   best_objective);
          break;
        }
        tree_.branch(node);
        break;
      case HighsMipStatus::kNodeInfeasible:
        reportMipSolverProgress(node_solve_status);
        if (options_.mip_report_level > 1) {
          printf("Node %9d (branch on %2d) infeasible\n", node.id,
                 node.branch_col);
          /*
            std::cout << "Node " << node.id
            << " infeasible." << std::endl;
          */
        }
        tree_.pop();
        break;
      case HighsMipStatus::kTimeout:
        return HighsMipStatus::kTimeout;
      case HighsMipStatus::kReachedSimplexIterationLimit:
        return HighsMipStatus::kReachedSimplexIterationLimit;
      case HighsMipStatus::kNodeUnbounded:
        return HighsMipStatus::kNodeUnbounded;
      default:
        /*
          std::cout << "Error or warning: Node " << node.id
          << " not solved to optimality, infeasibility or unboundedness." <<
          std::endl;
        */
        printf(
            "Node %9d (branch on %2d) not solved to optimality, infeasibility "
            "or unboundedness: status = %s\n",
            node.id, node.branch_col,
            highsMipStatusToString(node_solve_status).c_str());
        printf(
            "  Scaled model status is %s: max unscaled ( primal / dual ) "
            "infeasibilities are ( %g / %g )\n",
            highsModelStatusToString(scaled_model_status_).c_str(),
            info_.max_primal_infeasibility, info_.max_dual_infeasibility);
        printf("Unscaled model status is %s\n",
               highsModelStatusToString(model_status_).c_str());
        // Was break; but this causes infinite loop
        return HighsMipStatus::kNodeError;
    }
  }
  return HighsMipStatus::kTreeExhausted;
}

void HighsMipSolver::reportMipSolverProgress(const HighsMipStatus mip_status) {
  if (options_.mip_report_level == 1) {
    int report_frequency = 100;
    if (num_nodes_solved < 1000) {
      report_frequency = 100;
    } else if (num_nodes_solved < 10000) {
      report_frequency = 1000;
    } else if (num_nodes_solved < 100000) {
      report_frequency = 10000;
    } else {
      report_frequency = 100000;
    }
    switch (mip_status) {
      case HighsMipStatus::kOptimal:
        reportMipSolverProgressLine("");
        break;
      case HighsMipStatus::kTimeout:
        reportMipSolverProgressLine("Timeout");
        break;
      case HighsMipStatus::kReachedSimplexIterationLimit:
        reportMipSolverProgressLine("Reached simplex iteration limit");
        break;
      case HighsMipStatus::kError:
        reportMipSolverProgressLine("Error");
        break;
      case HighsMipStatus::kNodeOptimal:
        if (num_nodes_solved % report_frequency == 0)
          reportMipSolverProgressLine("");
        break;
      case HighsMipStatus::kNodeInfeasible:
        if (num_nodes_solved % report_frequency == 0)
          reportMipSolverProgressLine("");
        break;
      case HighsMipStatus::kNodeUnbounded:
        reportMipSolverProgressLine("Unbounded");
        break;
      case HighsMipStatus::kNodeNotOptimal:
        reportMipSolverProgressLine("Not optimal");
        break;
      case HighsMipStatus::kNodeError:
        reportMipSolverProgressLine("Node error");
        break;
      case HighsMipStatus::kRootNodeOptimal:
        reportMipSolverProgressLine("");
        break;
      case HighsMipStatus::kRootNodeNotOptimal:
        reportMipSolverProgressLine("Root node not optimal");
        break;
      case HighsMipStatus::kRootNodeError:
        reportMipSolverProgressLine("Root node error");
        break;
      case HighsMipStatus::kMaxNodeReached:
        reportMipSolverProgressLine("Max node reached");
        break;
      case HighsMipStatus::kUnderDevelopment:
        reportMipSolverProgressLine("Under development");
        break;
      case HighsMipStatus::kTreeExhausted:
        reportMipSolverProgressLine("Tree exhausted");
        break;
      default:
        reportMipSolverProgressLine("Unknown");
        break;
    }
  } else if (options_.mip_report_level > 1) {
    printf("Nodes solved = %d; Simplex iterations = %d\n", num_nodes_solved,
           info_.simplex_iteration_count);
  }
}

void HighsMipSolver::reportMipSolverProgressLine(std::string message,
                                                 const bool header) {
  if (header) {
    printf(
        "  Time |      Node |      Left |   LP iter | LP it/n |    dualbound | "
        " primalbound |    gap \n");
  } else {
    double average_simplex_iterations = info_.simplex_iteration_count;
    average_simplex_iterations /= num_nodes_solved;
    double time = timer_.readRunHighsClock();
    int num_nodes_left = tree_.getNumNodesLeft();
    double best_bound;
    double best_objective = tree_.getBestObjective();
    if (num_nodes_left > 0) {
      int best_node;
      best_bound = tree_.getBestBound(best_node);
    } else if (num_nodes_solved == 1) {
      // No nodes formed, so have just solved the root node
      best_bound = root_objective_;
      num_nodes_left = 2;
    } else {
      // No nodes left
      best_bound = best_objective;
      num_nodes_left = 0;
    }
    printf("%6.1f | %9d | %9d | %9d | %7.2f ", time, num_nodes_solved,
           num_nodes_left, info_.simplex_iteration_count,
           average_simplex_iterations);
    if (best_bound < HIGHS_CONST_INF) {
      printf("| %12.5e ", best_bound);
    } else {
      printf("|      --      ");
    }
    if (best_objective < HIGHS_CONST_INF) {
      // There is an integer solution
      double gap =
          100 * (best_objective - best_bound) / max(1.0, fabs(best_objective));
      printf("| %12.5e | %6.2f%%", best_objective, gap);
    } else {
      printf("|      --      |    Inf ");
    }
    printf(" %s\n", message.c_str());
  }
}

std::string HighsMipSolver::highsMipStatusToString(
    const HighsMipStatus mip_status) {
  switch (mip_status) {
    case HighsMipStatus::kOptimal:
      return "Optimal";
      break;
    case HighsMipStatus::kTimeout:
      return "Timeout";
      break;
    case HighsMipStatus::kReachedSimplexIterationLimit:
      return "Reached simplex iteration limit";
      break;
    case HighsMipStatus::kError:
      return "Error";
      break;
    case HighsMipStatus::kNodeOptimal:
      return "Node optimal";
      break;
    case HighsMipStatus::kNodeInfeasible:
      return "Node infeasible";
      break;
    case HighsMipStatus::kNodeUnbounded:
      return "Node unbounded";
      break;
    case HighsMipStatus::kNodeNotOptimal:
      return "Node not optimal";
      break;
    case HighsMipStatus::kNodeError:
      return "Node error";
      break;
    case HighsMipStatus::kRootNodeNotOptimal:
      return "Root node not optimal";
      break;
    case HighsMipStatus::kRootNodeError:
      return "Root node error";
      break;
    case HighsMipStatus::kMaxNodeReached:
      return "Max node number reached";
      break;
    case HighsMipStatus::kUnderDevelopment:
      return "Under development";
      break;
    case HighsMipStatus::kTreeExhausted:
      return "Tree exhausted";
      break;
    default:
#ifdef HiGHSDEV
      printf("HiGHS MIP status %d not recognised\n", (int)mip_status);
#endif
      return "Unrecognised HiGHS MIP status";
      break;
  }
  return "";
}
