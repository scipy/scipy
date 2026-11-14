/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_QPSOLVER_ASM_HPP__
#define __SRC_LIB_QPSOLVER_ASM_HPP__

#include "lp_data/HighsCallback.h"
#include "qpsolver/instance.hpp"
#include "qpsolver/qpconst.hpp"
#include "qpsolver/settings.hpp"
#include "qpsolver/statistics.hpp"
#include "util/HighsTimer.h"

enum class QpAsmStatus {
  kOk,
  kWarning,
  kError
  //  NEGATIVEEIGENVALUEINREDUCEDHESSIAN,
  //  BASISRANKDEFICIENT
};

struct QpSolution {
  QpVector primal;
  QpVector rowactivity;
  QpVector dualvar;
  QpVector dualcon;

  std::vector<BasisStatus> status_var;
  std::vector<BasisStatus> status_con;

  QpSolution(Instance& instance)
      : primal(QpVector(instance.num_var)),
        rowactivity(QpVector(instance.num_con)),
        dualvar(instance.num_var),
        dualcon(instance.num_con),
        status_var(instance.num_var),
        status_con(instance.num_con) {}
};

struct QpHotstartInformation {
  std::vector<HighsInt> active;
  std::vector<HighsInt> inactive;
  std::vector<BasisStatus> status;
  QpVector primal;
  QpVector rowact;

  QpHotstartInformation(HighsInt num_var, HighsInt num_row)
      : primal(QpVector(num_var)), rowact(QpVector(num_row)) {}
};

// the purpose of this is the pure algorithmic solution of a QP instance with
// given hotstart information. scenarios: 1) start from a given phase1 solution
// 2) start from a user-given hotstart solution
// 3) start from a qp solution that was attained from a scaled instance and
// cleanup 4) start from a qp solution that was attained from a perturbed
// instance and cleanup 5) start from a qp solution and cleanup after
// recomputing basis and reduced hessian factorization

std::string qpBasisStatusToString(const BasisStatus qp_basis_status);
std::string qpModelStatusToString(const QpModelStatus qp_model_status);
void assessQpPrimalFeasibility(
    const Instance& instance, const double primal_feasibility_tolerance,
    const std::vector<double>& var_value, const std::vector<double>& con_value,
    HighsInt& num_var_infeasibilities, double& max_var_infeasibility,
    double& sum_var_infeasibilities, HighsInt& num_con_infeasibilities,
    double& max_con_infeasibility, double& sum_con_infeasibilities,
    double& max_con_residual, double& sum_con_residuals);

QpAsmStatus solveqpActual(Instance& instance, Settings& settings,
                          QpHotstartInformation& startinfo, Statistics& stats,
                          QpModelStatus& status, QpSolution& solution,
                          HighsTimer& qp_timer, HighsCallback& callback);

#endif
