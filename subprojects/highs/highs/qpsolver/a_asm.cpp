/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "qpsolver/a_asm.hpp"

#include "qpsolver/quass.hpp"
#include "util/HighsCDouble.h"

QpAsmStatus solveqpActual(Instance& instance, Settings& settings,
                          QpHotstartInformation& startinfo, Statistics& stats,
                          QpModelStatus& status, QpSolution& solution,
                          HighsTimer& qp_timer, HighsCallback& callback) {
  Runtime rt(instance, stats);
  rt.settings = settings;
  Quass quass(rt);

  Basis basis(rt, startinfo.active, startinfo.status, startinfo.inactive);

  quass.solve(startinfo.primal, startinfo.rowact, basis, qp_timer, callback);

  status = rt.status;

  solution.status_var = rt.status_var;
  solution.status_con = rt.status_con;
  solution.primal = rt.primal;
  solution.dualvar = rt.dualvar;
  solution.rowactivity = rt.rowactivity;
  solution.dualcon = rt.dualcon;

  return QpAsmStatus::kOk;
}

std::string qpBasisStatusToString(const BasisStatus qp_basis_status) {
  switch (qp_basis_status) {
    case BasisStatus::kInactive:
      return "Inactive";
    case BasisStatus::kActiveAtLower:
      return "Active at lower bound";
    case BasisStatus::kActiveAtUpper:
      return "Active at upper bound";
    case BasisStatus::kInactiveInBasis:
      return "Inactive in basis";
    default:
      return "Unidentified QP basis status";
  }
}

std::string qpModelStatusToString(const QpModelStatus qp_model_status) {
  switch (qp_model_status) {
    case QpModelStatus::kNotset:
      return "Not set";
    case QpModelStatus::kUndetermined:
      return "Undetermined";
    case QpModelStatus::kOptimal:
      return "Optimal";
    case QpModelStatus::kUnbounded:
      return "Unbounded";
    case QpModelStatus::kInfeasible:
      return "Infeasible";
    case QpModelStatus::kIterationLimit:
      return "Iteration limit";
    case QpModelStatus::kTimeLimit:
      return "Time ;limit";
    case QpModelStatus::kLargeNullspace:
      return "Large nullspace";
    case QpModelStatus::kNonConvex:
      return "Non-convex";
    case QpModelStatus::kError:
      return "Error";
    default:
      return "Unidentified QP model status";
  }
}

void assessQpPrimalFeasibility(
    const Instance& instance, const double primal_feasibility_tolerance,
    const std::vector<double>& var_value, const std::vector<double>& con_value,
    HighsInt& num_var_infeasibilities, double& max_var_infeasibility,
    double& sum_var_infeasibilities, HighsInt& num_con_infeasibilities,
    double& max_con_infeasibility, double& sum_con_infeasibilities,
    double& max_con_residual, double& sum_con_residuals) {
  num_var_infeasibilities = 0;
  max_var_infeasibility = 0;
  sum_var_infeasibilities = 0;
  num_con_infeasibilities = 0;
  max_con_infeasibility = 0;
  sum_con_infeasibilities = 0;
  max_con_residual = 0;
  sum_con_residuals = 0;
  // Valid solution, but is it feasible?
  std::vector<HighsCDouble> con_value_quad;
  con_value_quad.assign(instance.num_con, HighsCDouble{0.0});
  for (HighsInt iVar = 0; iVar < instance.num_var; iVar++) {
    double lower = instance.var_lo[iVar];
    double upper = instance.var_up[iVar];
    double primal = var_value[iVar];
    double var_infeasibility = 0;
    if (primal < lower - primal_feasibility_tolerance) {
      var_infeasibility = lower - primal;
    } else if (primal > upper + primal_feasibility_tolerance) {
      var_infeasibility = primal - upper;
    }
    if (var_infeasibility > 0) {
      if (var_infeasibility > primal_feasibility_tolerance)
        num_var_infeasibilities++;
      max_var_infeasibility =
          std::max(var_infeasibility, max_var_infeasibility);
      sum_var_infeasibilities += var_infeasibility;
    }
    for (HighsInt iEl = instance.A.mat.start[iVar];
         iEl < instance.A.mat.start[iVar + 1]; iEl++) {
      con_value_quad[instance.A.mat.index[iEl]] +=
          primal * instance.A.mat.value[iEl];
    }
  }
  for (HighsInt iCon = 0; iCon < instance.num_con; iCon++) {
    double lower = instance.con_lo[iCon];
    double upper = instance.con_up[iCon];
    double primal = con_value[iCon];
    double con_infeasibility = 0;
    if (primal < lower - primal_feasibility_tolerance) {
      con_infeasibility = lower - primal;
    } else if (primal > upper + primal_feasibility_tolerance) {
      con_infeasibility = primal - upper;
    }
    if (con_infeasibility > 0) {
      if (con_infeasibility > primal_feasibility_tolerance)
        num_con_infeasibilities++;
      max_con_infeasibility =
          std::max(con_infeasibility, max_con_infeasibility);
      sum_con_infeasibilities += con_infeasibility;
    }
    double con_residual = std::fabs(primal - double(con_value_quad[iCon]));
    max_con_residual = std::max(con_residual, max_con_residual);
    sum_con_residuals += con_residual;
  }
}
