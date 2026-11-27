/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file qpsolver/QpAsmWrapper.cpp
 * @brief
 */

#include "qpsolver/QpAsmWrapper.h"

#include "model/HighsHessianUtils.h"
#include "qpsolver/a_quass.hpp"
#include "qpsolver/runtime.hpp"

HighsStatus solveQpAsm(HighsQpSolverObject& solver_object) {
  return solveQpAsm(solver_object.options_, solver_object.timer_,
                    solver_object.model_.lp_, solver_object.model_.hessian_,
                    solver_object.basis_, solver_object.solution_,
                    solver_object.model_status_, solver_object.highs_info_,
                    solver_object.callback_);
}

HighsHessianFunctionType testOracleCallSquareHessian =
    [](const HighsInt call_type, const HighsInt* x_num_entries,
       const HighsInt* x_index, const double* x_value,
       HighsInt* hessian_x_num_entries, HighsInt* hessian_x_index,
       double* hessian_x_value, void* hessian_p) {
      assert(kHessianOracleCallTypeMin <= call_type &&
             call_type <= kHessianOracleCallTypeMax);

      HighsHessian hessian = *(static_cast<HighsHessian*>(hessian_p));
      assert(hessian.format_ == HessianFormat::kSquare);

      // Lambda for adding multiple of Hessian column into hessian_x_value
      auto addScaledQcol = [&](const HighsInt iCol, const double x_value) {
        for (HighsInt iEl = hessian.start_[iCol];
             iEl < hessian.start_[iCol + 1]; iEl++) {
          HighsInt iRow = hessian.index_[iEl];
          hessian_x_value[iRow] += hessian.value_[iEl] * x_value;
        }
      };

      if (call_type == kHessianOracleCallTypeEntry) {
        assert(x_num_entries == nullptr);
        assert(x_value == nullptr);
        assert(x_index != nullptr);
        assert(hessian_x_num_entries == nullptr);
        assert(hessian_x_index != nullptr);
        assert(hessian_x_value != nullptr);
        HighsInt iCol = x_index[0];
        HighsInt iRow = hessian_x_index[0];
        // Zero Qx value in case the Hessian entry requested is zero
        hessian_x_value[0] = 0;
        for (HighsInt iEl = hessian.start_[iCol];
             iEl < hessian.start_[iCol + 1]; iEl++) {
          if (hessian.index_[iEl] == iRow) {
            hessian_x_value[0] = hessian.value_[iEl];
            return 0;
          }
        }
      } else if (call_type == kHessianOracleCallTypeColumn) {
        // Get the entries in column iCol
        assert(x_num_entries == nullptr);
        assert(x_value == nullptr);
        assert(x_index != nullptr);
        assert(hessian_x_num_entries != nullptr);
        assert(hessian_x_index != nullptr);
        assert(hessian_x_value != nullptr);
        (*hessian_x_num_entries) = 0;
        HighsInt iCol = x_index[0];
        for (HighsInt iEl = hessian.start_[iCol];
             iEl < hessian.start_[iCol + 1]; iEl++) {
          hessian_x_index[*hessian_x_num_entries] = hessian.index_[iEl];
          hessian_x_value[*hessian_x_num_entries] = hessian.value_[iEl];
          (*hessian_x_num_entries)++;
        }
      } else {
        assert(x_index == nullptr || *x_num_entries >= 0);
        assert(hessian_x_num_entries == nullptr);
        assert(hessian_x_index == nullptr);
        assert(hessian_x_value != nullptr);
        if (x_index == nullptr) {
          // Simple product with full vector x, full vector q_x
          for (HighsInt iCol = 0; iCol < hessian.dim_; iCol++)
            addScaledQcol(iCol, x_value[iCol]);
        } else {
          // x is scattered with x_num_entries entries in rows x_index
          for (HighsInt iX = 0; iX < *x_num_entries; iX++) {
            HighsInt iCol = x_index[iX];
            addScaledQcol(iCol, x_value[iCol]);
          }
        }
      }
      return 0;
    };

HighsStatus solveQpAsm(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, const HighsHessian& hessian,
                       HighsBasis& basis, HighsSolution& solution,
                       HighsModelStatus& model_status, HighsInfo& info,
                       HighsCallback& callback) {
  Instance instance(lp.num_col_, lp.num_row_);

  instance.sense = HighsInt(lp.sense_);
  instance.num_con = lp.num_row_;
  instance.num_var = lp.num_col_;

  instance.A.mat.num_col = lp.num_col_;
  instance.A.mat.num_row = lp.num_row_;
  instance.A.mat.start = lp.a_matrix_.start_;
  instance.A.mat.index = lp.a_matrix_.index_;
  instance.A.mat.value = lp.a_matrix_.value_;
  instance.c.value = lp.col_cost_;
  instance.offset = lp.offset_;
  instance.con_lo = lp.row_lower_;
  instance.con_up = lp.row_upper_;
  instance.var_lo = lp.col_lower_;
  instance.var_up = lp.col_upper_;
  // Clear the instance Hessian data
  instance.Q.mat.clear();
  HighsHessian oracle_hessian;
  if (hessian.dim_ > 0) {
    if (options.test_qp_oracle) {
      // Test the Hessian oracle by using the incumbent Hessian as data for it
      oracle_hessian = hessian.toSquare();
      instance.Q.mat.oracle_.dim_ = hessian.dim_;
      instance.Q.mat.oracle_.call_ = testOracleCallSquareHessian;
      instance.Q.mat.oracle_.data_ = &oracle_hessian;
    } else {
      assert(instance.Q.mat.oracle_.call_ == nullptr);
    }
    // Generate the explicit Hessian data to use if not testing the
    // oracle, and to cross-check the oracle results
    instance.Q.mat.num_col = lp.num_col_;
    instance.Q.mat.num_row = lp.num_col_;
    triangularToSquareHessian(hessian, instance.Q.mat.start,
                              instance.Q.mat.index, instance.Q.mat.value);
  } else {
    assert(hessian.isOracle());
    instance.Q.mat.oracle_ = hessian.oracle_;
  }

  for (HighsInt i = 0; i < (HighsInt)instance.c.value.size(); i++) {
    if (instance.c.value[i] != 0.0) {
      instance.c.index[instance.c.num_nz++] = i;
    }
  }

  if (lp.sense_ == ObjSense::kMaximize) {
    // Negate the vector and Hessian
    for (double& i : instance.c.value) i *= -1.0;
    if (instance.Q.mat.num_col > 0) {
      for (double& i : instance.Q.mat.value) i *= -1.0;
    }
    if (instance.Q.mat.isOracle()) instance.Q.mat.oracle_.multiplier_ = -1;
  }

  Settings settings;
  Statistics stats;

  settings.reportingfequency = 100;
  if (options.qp_iteration_limit <= 10) {
    settings.reportingfequency = 1;
  } else if (options.qp_iteration_limit <= 100) {
    settings.reportingfequency = 10;
  }
  // Setting qp_update_limit = 10 leads to error with lpHighs3
  const HighsInt qp_update_limit = 1000;  // 1000; // default
  if (qp_update_limit != settings.reinvertfrequency) {
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "Changing QP reinversion frequency from %d to %d\n",
                 int(settings.reinvertfrequency), int(qp_update_limit));
    settings.reinvertfrequency = qp_update_limit;
  }

  settings.allow_hot_start = options.qp_allow_hot_start;
  settings.iteration_limit = options.qp_iteration_limit;
  settings.nullspace_limit = options.qp_nullspace_limit;
  assert(settings.hessian_regularization_value == kHessianRegularizationValue);
  settings.hessian_regularization_value = options.qp_regularization_value;
  settings.primal_feasibility_tolerance = options.primal_feasibility_tolerance;

  // Define the QP model status logging function
  settings.qp_model_status_log.subscribe([&](QpModelStatus& qp_model_status) {
    if (qp_model_status == QpModelStatus::kUndetermined ||
        qp_model_status == QpModelStatus::kLargeNullspace ||
        qp_model_status == QpModelStatus::kNonConvex ||
        qp_model_status == QpModelStatus::kError ||
        qp_model_status == QpModelStatus::kNotset)
      highsLogUser(options.log_options, HighsLogType::kInfo,
                   "QP solver model status: %s\n",
                   qpModelStatusToString(qp_model_status).c_str());
  });

  // Define the QP solver iteration logging function
  settings.iteration_log_header.subscribe([&](HighsInt& null) {
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "  Iteration        Objective     NullspaceDim\n");
  });

  // Define the QP solver iteration logging function
  settings.iteration_log.subscribe([&](Statistics& stats) {
    int rep = stats.iteration.size() - 1;
    std::string time_string =
        options.timeless_log ? ""
                             : highsFormatToString(" %9.2fs", stats.time[rep]);
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "%11d  %15.8g           %6d%s\n", int(stats.iteration[rep]),
                 stats.objval[rep], int(stats.nullspacedimension[rep]),
                 time_string.c_str());
  });

  // Define the QP nullspace limit logging function
  settings.nullspace_limit_log.subscribe([&](HighsInt& nullspace_limit) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "QP solver has exceeded nullspace limit of %d\n",
                 int(nullspace_limit));
  });

  // Define the degeneracy failure logging function
  settings.degeneracy_fail_log.subscribe(
      [&](std::pair<HighsInt, double>& degeneracy_fail_data) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "QP solver has failed due to degeneracy: "
                     "cannot find non-active constraint to leave basis."
                     " max: log(d[%d]) = %lf\n",
                     int(degeneracy_fail_data.first),
                     degeneracy_fail_data.second);
      });

  settings.time_limit = options.time_limit;
  settings.lambda_zero_threshold = options.dual_feasibility_tolerance;

  switch (options.simplex_primal_edge_weight_strategy) {
    case 0:
      settings.pricing = PricingStrategy::DantzigWolfe;
      break;
    case 1:
      settings.pricing = PricingStrategy::Devex;
      break;
    case 2:
      settings.pricing = PricingStrategy::SteepestEdge;
      break;
    default:
      settings.pricing = PricingStrategy::Devex;
  }

  QpAsmStatus status = solveqp(instance, settings, stats, model_status, basis,
                               solution, timer, callback);

  // QP solver can fail, so should return something other than
  // QpAsmStatus::kOk
  if (status == QpAsmStatus::kError) return HighsStatus::kError;

  assert(status == QpAsmStatus::kOk || status == QpAsmStatus::kWarning);
  HighsStatus return_status = status == QpAsmStatus::kWarning
                                  ? HighsStatus::kWarning
                                  : HighsStatus::kOk;

  // Set the QP-specific values of info
  info.simplex_iteration_count += stats.phase1_iterations;
  info.qp_iteration_count += stats.num_iterations;
  // Basis validity is set locally to solvers; solution status is set
  // in getKktErrors
  info.basis_validity =
      basis.valid ? kBasisValidityValid : kBasisValidityInvalid;
  return return_status;
}
