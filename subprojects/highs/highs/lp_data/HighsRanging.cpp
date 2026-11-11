/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsRanging.cpp
 * @brief Compute LP ranging data for HiGHS
 */
#include "lp_data/HighsRanging.h"

#include <algorithm>
#include <cassert>
#include <functional>  // for negate
#include <sstream>

#include "lp_data/HighsModelUtils.h"

using std::min;

void HighsRanging::invalidate() { valid = false; }

void HighsRanging::clear() {
  this->invalidate();
  this->col_cost_up.value_.clear();
  this->col_cost_up.objective_.clear();
  this->col_cost_up.in_var_.clear();
  this->col_cost_up.ou_var_.clear();
  this->col_cost_dn.value_.clear();
  this->col_cost_dn.objective_.clear();
  this->col_cost_dn.in_var_.clear();
  this->col_cost_dn.ou_var_.clear();
  this->col_bound_up.value_.clear();
  this->col_bound_up.objective_.clear();
  this->col_bound_up.in_var_.clear();
  this->col_bound_up.ou_var_.clear();
  this->col_bound_dn.value_.clear();
  this->col_bound_dn.objective_.clear();
  this->col_bound_dn.in_var_.clear();
  this->col_bound_dn.ou_var_.clear();
  this->row_bound_up.value_.clear();
  this->row_bound_up.objective_.clear();
  this->row_bound_up.in_var_.clear();
  this->row_bound_up.ou_var_.clear();
  this->row_bound_dn.value_.clear();
  this->row_bound_dn.objective_.clear();
  this->row_bound_dn.in_var_.clear();
  this->row_bound_dn.ou_var_.clear();
}

static double infProduct(double value) {
  // Multiplying value and kHighsInf
  if (value == 0) {
    return 0;
  } else {
    return value * kHighsInf;
  }
}

static double possInfProduct(double poss_inf, double value) {
  // Multiplying something that could be infinite and value
  if (value == 0) {
    return 0;
  } else {
    return poss_inf * value;
  }
}

HighsStatus getRangingData(HighsRanging& ranging,
                           HighsLpSolverObject& solver_object) {
  ranging.clear();
  if (solver_object.model_status_ != HighsModelStatus::kOptimal) {
    highsLogUser(solver_object.options_.log_options, HighsLogType::kError,
                 "Cannot get ranging without an optimal solution\n");
    return HighsStatus::kError;
  }
  HEkk& ekk_instance = solver_object.ekk_instance_;
  if (!ekk_instance.status_.initialised_for_solve) {
    highsLogUser(solver_object.options_.log_options, HighsLogType::kError,
                 "Cannot get ranging without a valid Simplex instance\n");
    return HighsStatus::kError;
  }
  const HighsLp& use_lp = solver_object.lp_;

  // Unscale the simplex data if the LP has been solved in the scaled space
  ekk_instance.unscaleSimplex(use_lp);

  // Aliases
  const HighsSimplexInfo& simplex_info = ekk_instance.info_;
  const SimplexBasis& simplex_basis = ekk_instance.basis_;
  const vector<double>& value_ = simplex_info.workValue_;
  const vector<double>& dual_ = simplex_info.workDual_;
  const vector<double>& cost_ = simplex_info.workCost_;
  const vector<double>& lower_ = simplex_info.workLower_;
  const vector<double>& upper_ = simplex_info.workUpper_;
  const vector<double>& Bvalue_ = simplex_info.baseValue_;
  const vector<double>& Blower_ = simplex_info.baseLower_;
  const vector<double>& Bupper_ = simplex_info.baseUpper_;
  const vector<int8_t>& Nflag_ = simplex_basis.nonbasicFlag_;
  const vector<int8_t>& Nmove_ = simplex_basis.nonbasicMove_;
  const vector<HighsInt>& Bindex_ = simplex_basis.basicIndex_;
  const HighsSparseMatrix& matrix = use_lp.a_matrix_;

  // Local copies of scalars

  const HighsInt numRow = use_lp.num_row_;
  const HighsInt numCol = use_lp.num_col_;
  const HighsInt numTotal = numCol + numRow;
  const double H_TT = 1e-13;
  const double H_INF = kHighsInf;
  const double objective = solver_object.highs_info_.objective_function_value;

  // Code written for minimization problems. Maximization problems are
  // solved by using negated costs in the simplex solver and
  // minimizing. Thus dual information in the simplex solver is
  // negated for maximization problems. The objective has the right
  // sign, though. Maximization problems are, thus, accommodated by
  // applying the sign multiplier to dual information.
  HighsInt sense = 1;
  if (use_lp.sense_ == ObjSense::kMaximize) sense = -1;

  vector<HighsInt> iWork_(numTotal);
  vector<double> dWork_(numTotal);
  HVector column;
  column.setup(numRow);

  vector<double> xi = Bvalue_;
  for (HighsInt i = 0; i < numRow; i++) {
    xi[i] = max(xi[i], Blower_[i]);
    xi[i] = min(xi[i], Bupper_[i]);
  }
  vector<double> dj = dual_;
  for (HighsInt j = 0; j < numTotal; j++) {
    if (Nflag_[j] && (lower_[j] != upper_[j])) {
      if (value_[j] == lower_[j]) dj[j] = max(dj[j], 0.0);
      if (value_[j] == upper_[j]) dj[j] = min(dj[j], 0.0);
      if (lower_[j] == -H_INF && upper_[j] == H_INF) dj[j] = 0;
    }
  }
  //
  // Ranging 1.2. prepare "delta" space
  //
  vector<double> dxi_inc(numRow);
  vector<double> dxi_dec(numRow);
  for (HighsInt i = 0; i < numRow; i++) {
    dxi_inc[i] = Bupper_[i] - xi[i];
    dxi_dec[i] = Blower_[i] - xi[i];
  }

  vector<double> ddj_inc(numTotal);
  vector<double> ddj_dec(numTotal);
  for (HighsInt j = 0; j < numTotal; j++) {
    if (Nflag_[j]) {
      ddj_inc[j] = (value_[j] == lower_[j]) ? +H_INF : -dj[j];
      ddj_dec[j] = (value_[j] == upper_[j]) ? -H_INF : -dj[j];
    }
  }

  //
  // Ranging 1.3. prepare "theta" space
  //
  const double tol_a = 1e-9;
  const double THETA_INF = H_INF / 1e40;

  vector<double> txj_inc(numTotal, +THETA_INF);  // theta
  vector<double> axj_inc(numTotal, 0);           // alpha
  vector<HighsInt> ixj_inc(numTotal, -1);        // i-out
  vector<HighsInt> wxj_inc(numTotal, 0);         // which bound is limiting
  vector<HighsInt> jxj_inc(numTotal, -1);        // j = n(i), (with bound flip)

  vector<double> txj_dec(numTotal, -THETA_INF);
  vector<double> axj_dec(numTotal, 0);
  vector<HighsInt> ixj_dec(numTotal, -1);
  vector<HighsInt> wxj_dec(numTotal, 0);
  vector<HighsInt> jxj_dec(numTotal, -1);

  vector<double> tci_inc(numRow, +THETA_INF);  // theta
  vector<double> aci_inc(numRow, 0);           // alpha
  vector<HighsInt> jci_inc(numRow, -1);        // column index

  vector<double> tci_dec(numRow, -THETA_INF);
  vector<double> aci_dec(numRow, 0);
  vector<HighsInt> jci_dec(numRow, -1);

  // Major "theta" loop
  for (HighsInt j = 0; j < numTotal; j++) {
    // Skip basic column
    if (!Nflag_[j]) continue;

    // Form updated column
    column.clear();
    matrix.collectAj(column, j, 1);

    const double expected_density = ekk_instance.info_.col_aq_density;
    ekk_instance.ftran(column, expected_density);
    HighsInt nWork = 0;
    for (HighsInt k = 0; k < column.count; k++) {
      HighsInt iRow = column.index[k];
      double alpha = column.array[iRow];
      if (fabs(alpha) > tol_a) {
        iWork_[nWork] = iRow;
        dWork_[nWork] = alpha;
        nWork++;
      }
    }
    // Standard primal ratio test
    double myt_inc = +THETA_INF;
    double myt_dec = -THETA_INF;
    HighsInt myk_inc = -1;
    HighsInt myk_dec = -1;
    for (HighsInt k = 0; k < nWork; k++) {
      HighsInt i = iWork_[k];
      double alpha = dWork_[k];
      double theta_inc = (alpha < 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      double theta_dec = (alpha > 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      if (myt_inc > theta_inc) myt_inc = theta_inc, myk_inc = k;
      if (myt_dec < theta_dec) myt_dec = theta_dec, myk_dec = k;
    }

    if (myk_inc != -1) {
      HighsInt i = iWork_[myk_inc];
      double alpha = dWork_[myk_inc];
      ixj_inc[j] = i;
      axj_inc[j] = alpha;
      txj_inc[j] = (alpha < 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      wxj_inc[j] = (alpha < 0 ? +1 : -1);
    }

    if (myk_dec != -1) {
      HighsInt i = iWork_[myk_dec];
      double alpha = dWork_[myk_dec];
      ixj_dec[j] = i;
      axj_dec[j] = alpha;
      txj_dec[j] = (alpha > 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      wxj_dec[j] = (alpha > 0 ? +1 : -1);
    }

    // Accumulated dual ratio test
    double myd_inc = ddj_inc[j];
    double myd_dec = ddj_dec[j];
    for (HighsInt k = 0; k < nWork; k++) {
      HighsInt i = iWork_[k];
      double alpha = dWork_[k];
      double theta_inc = (alpha < 0 ? myd_inc : myd_dec) / -alpha;
      double theta_dec = (alpha > 0 ? myd_inc : myd_dec) / -alpha;
      if (tci_inc[i] > theta_inc)
        tci_inc[i] = theta_inc, aci_inc[i] = alpha, jci_inc[i] = j;
      if (tci_dec[i] < theta_dec)
        tci_dec[i] = theta_dec, aci_dec[i] = alpha, jci_dec[i] = j;
    }
  }

  // Additional j-out for primal ratio test (considering bound flip)
  for (HighsInt j = 0; j < numTotal; j++) {
    if (Nflag_[j]) {
      // J-out for x_j = l_j
      if (Nmove_[j] == +1) {
        double value = value_[j] + txj_inc[j];
        if (ixj_inc[j] != -1 && value <= upper_[j]) {
          jxj_inc[j] = Bindex_[ixj_inc[j]];
        } else if (value > upper_[j]) {
          jxj_inc[j] = j;
        }
      }
      // J-out for x_j = u_j
      if (Nmove_[j] == -1) {
        double value = value_[j] + txj_dec[j];
        if (ixj_dec[j] != -1 && value >= lower_[j]) {
          jxj_dec[j] = Bindex_[ixj_dec[j]];
        } else if (value < lower_[j]) {
          jxj_dec[j] = j;
        }
      }
      // J-out for free variable
      if (lower_[j] == -H_INF && upper_[j] == H_INF) {
        if (ixj_inc[j] != -1) jxj_inc[j] = jxj_dec[j] = Bindex_[ixj_inc[j]];
        if (ixj_dec[j] != -1) jxj_inc[j] = jxj_dec[j] = Bindex_[ixj_dec[j]];
      }
    }
  }

  //
  // Ranging 2. cost ranging
  //
  vector<double> c_up_c(numTotal), c_dn_c(numTotal);
  vector<double> c_up_f(numTotal), c_dn_f(numTotal);
  vector<HighsInt> c_up_e(numTotal), c_dn_e(numTotal);
  vector<HighsInt> c_up_l(numTotal), c_dn_l(numTotal);

  //
  // Ranging 2.1. nonbasic cost ranging
  //
  //  const HighsInt check_col = 2951;
  for (HighsInt j = 0; j < numCol; j++) {
    if (Nflag_[j]) {
      // Primal value and its sign
      double value = value_[j];
      double vsign = (value > 0) ? 1 : (value < 0 ? -1 : 0);

      // Increase c_j
      if (ddj_inc[j] != H_INF) {
        c_up_c[j] = cost_[j] + ddj_inc[j];
        c_up_f[j] =
            objective +
            sense * possInfProduct(ddj_inc[j], value);  // value * ddj_inc[j];
        c_up_e[j] = j;
        c_up_l[j] = jxj_dec[j];
      } else {
        c_up_c[j] = H_INF;
        c_up_f[j] = objective + sense * infProduct(vsign);  // vsign * H_INF;
        c_up_e[j] = -1;
        c_up_l[j] = -1;
      }

      // Decrease c_j
      if (ddj_dec[j] != H_INF) {
        c_dn_c[j] = cost_[j] + ddj_dec[j];
        c_dn_f[j] =
            objective +
            sense * possInfProduct(ddj_dec[j], value);  // value * ddj_dec[j];
        c_dn_e[j] = j;
        c_dn_l[j] = jxj_inc[j];
      } else {
        c_up_c[j] = -H_INF;
        c_up_f[j] = objective - sense * infProduct(vsign);  // vsign * H_INF;
        c_up_e[j] = -1;
        c_up_l[j] = -1;
      }
    }
  }

  //
  // Ranging 2.2. basic cost ranging
  //

  for (HighsInt i = 0; i < numRow; i++) {
    if (Bindex_[i] < numCol) {
      // Primal variable and its sign
      HighsInt j = Bindex_[i], je;
      double value = xi[i];
      double vsign = (value > 0) ? 1 : (value < 0 ? -1 : 0);

      // Increase c_i
      if (jci_inc[i] != -1) {
        c_up_c[j] = cost_[j] + tci_inc[i];
        c_up_f[j] =
            objective +
            sense * possInfProduct(tci_inc[i], value);  // value * tci_inc[i];
        c_up_e[j] = je = jci_inc[i];
        c_up_l[j] = Nmove_[je] > 0 ? jxj_inc[je] : jxj_dec[je];
      } else {
        c_up_c[j] = H_INF;
        c_up_f[j] = objective + sense * infProduct(vsign);  // vsign * H_INF;
        c_up_e[j] = -1;
        c_up_l[j] = -1;
      }

      // Decrease c_i
      if (jci_dec[i] != -1) {
        c_dn_c[j] = cost_[j] + tci_dec[i];
        c_dn_f[j] =
            objective +
            sense * possInfProduct(tci_dec[i], value);  // value * tci_dec[i];
        c_dn_e[j] = je = jci_dec[i];
        c_dn_l[j] = Nmove_[je] > 0 ? jxj_inc[je] : jxj_dec[je];
      } else {
        c_dn_c[j] = -H_INF;
        c_dn_f[j] = objective - sense * infProduct(vsign);  // H_INF * vsign;
        c_dn_e[j] = -1;
        c_dn_l[j] = -1;
      }
    }
  }

  //
  // Ranging 3. bounds ranging
  //
  vector<double> b_up_b(numTotal), b_dn_b(numTotal);
  vector<double> b_up_f(numTotal), b_dn_f(numTotal);
  vector<HighsInt> b_up_e(numTotal), b_dn_e(numTotal);
  vector<HighsInt> b_up_l(numTotal), b_dn_l(numTotal);

  //
  // Ranging 3.1. nonbasic bounds ranging
  //
  for (HighsInt j = 0; j < numTotal; j++) {
    if (Nflag_[j]) {
      // FREE variable
      if (lower_[j] == -H_INF && upper_[j] == H_INF) {
        b_up_b[j] = H_INF;
        b_up_f[j] = objective;
        b_up_e[j] = -1;
        b_up_l[j] = -1;
        b_dn_b[j] = -H_INF;
        b_dn_f[j] = objective;
        b_dn_e[j] = -1;
        b_dn_l[j] = -1;
        continue;
      }

      // Dual value and its sign
      double dualv = dj[j];
      double dsign = (dualv > 0) ? 1 : (dualv < 0 ? -1 : 0);

      // Increase x_j
      if (ixj_inc[j] != -1) {
        HighsInt i = ixj_inc[j];
        b_up_b[j] = value_[j] + txj_inc[j];
        b_up_f[j] =
            objective +
            sense * possInfProduct(txj_inc[j], dualv);  // txj_inc[j] * dualv;
        b_up_e[j] = wxj_inc[j] > 0 ? jci_inc[i] : jci_dec[i];
        b_up_l[j] = Bindex_[i];
      } else {
        b_up_b[j] = H_INF;
        b_up_f[j] = objective + sense * infProduct(dsign);  // H_INF * dsign;
        b_up_e[j] = -1;
        b_up_l[j] = -1;
      }

      // Check if b_up_b > upper
      if (value_[j] != upper_[j] && b_up_b[j] > upper_[j]) {
        b_up_b[j] = upper_[j];
        assert(lower_[j] > -kHighsInf);
        b_up_f[j] = objective + sense * (upper_[j] - lower_[j]) * dualv;
        b_up_e[j] = j;
        b_up_l[j] = j;
      }

      // Decrease x_j
      if (ixj_dec[j] != -1) {
        HighsInt i = ixj_dec[j];
        b_dn_b[j] = value_[j] + txj_dec[j];
        b_dn_f[j] =
            objective +
            sense * possInfProduct(txj_dec[j], dualv);  // txj_dec[j] * dualv;
        b_dn_e[j] = wxj_dec[j] > 0 ? jci_inc[i] : jci_dec[i];
        b_dn_l[j] = Bindex_[i];
      } else {
        b_dn_b[j] = -H_INF;
        b_dn_f[j] = objective - sense * infProduct(dsign);  // H_INF * dsign;
        b_dn_e[j] = -1;
        b_dn_l[j] = -1;
      }

      // Check if b_dn_b < lower
      if (value_[j] != lower_[j] && b_dn_b[j] < lower_[j]) {
        b_dn_b[j] = lower_[j];
        assert(upper_[j] < kHighsInf);
        b_dn_f[j] = objective + sense * (lower_[j] - upper_[j]) * dualv;
        b_dn_e[j] = j;
        b_dn_l[j] = j;
      }
    }
  }

  //
  // Ranging 3.2. basic bounds ranging
  //
  for (HighsInt i = 0; i < numRow; i++) {
    for (HighsInt dir = -1; dir <= 1; dir += 2) {
      HighsInt j = Bindex_[i];
      double& newx = dir == -1 ? b_dn_b[j] : b_up_b[j];
      double& newf = dir == -1 ? b_dn_f[j] : b_up_f[j];
      HighsInt& j_enter = dir == -1 ? b_dn_e[j] : b_up_e[j];
      HighsInt& j_leave = dir == -1 ? b_dn_l[j] : b_up_l[j];

      HighsInt j_in = dir == -1 ? jci_inc[i] : jci_dec[i];
      double a_in = dir == -1 ? aci_inc[i] : aci_dec[i];
      if (j_in != -1) {
        HighsInt jmove = Nmove_[j_in];
        HighsInt i_out = jmove > 0 ? ixj_inc[j_in] : ixj_dec[j_in];
        HighsInt j_out = jmove > 0 ? jxj_inc[j_in] : jxj_dec[j_in];
        HighsInt w_out = jmove > 0 ? wxj_inc[j_in] : wxj_dec[j_in];
        double tt = jmove > 0 ? txj_inc[j_in] : txj_dec[j_in];
        if (j_out == j_in) {
          // Bound flip
          double delta = jmove * (upper_[j_in] - lower_[j_in]);
          newx = xi[i] - delta * a_in;
          newf = objective + sense * delta * dual_[j_in];
          j_enter = j_in;
          j_leave = j_out;
        } else if (j_out != -1) {
          // Regular
          double delta = w_out > 0 ? dxi_inc[i_out] : dxi_dec[i_out];
          double a_out = jmove > 0 ? axj_inc[j_in] : axj_dec[j_in];
          newx = xi[i] + delta * a_in / a_out;
          newf = objective + sense * tt * dual_[j_in];
          j_enter = j_in;
          j_leave = j_out;
        } else {
          // Primal ratio test failed - change unlimitedly
          //
          // While still limited by its own bounds
          //
          // Its own bounds could just be inf
          newx = dir == -1 ? lower_[j] : upper_[j];
          newf = objective;
          j_enter = -1;
          j_leave = -1;
        }
      } else {
        // Dual ratio test failed - just stay
        newx = xi[i];
        newf = objective;
        j_enter = -1;
        j_leave = -1;
      }
    }
  }

  //
  // Ranging 4.1.1 Trim small value to zero
  //

  for (HighsInt j = 0; j < numCol; j++) {
    if (fabs(c_up_c[j]) < H_TT) c_up_c[j] = 0;
    if (fabs(c_dn_c[j]) < H_TT) c_dn_c[j] = 0;
    if (fabs(b_up_b[j]) < H_TT) b_up_b[j] = 0;
    if (fabs(b_dn_b[j]) < H_TT) b_dn_b[j] = 0;
  }
  for (HighsInt i = 0, j = numCol; i < numRow; i++, j++) {
    if (fabs(b_up_b[j]) < H_TT) b_up_b[j] = 0;
    if (fabs(b_dn_b[j]) < H_TT) b_dn_b[j] = 0;
  }

  //
  // Ranging 4.2. Put to output buffer
  //

  if (sense > 0) {
    ranging.col_cost_up.value_ = c_up_c;
    ranging.col_cost_dn.value_ = c_dn_c;
    ranging.col_cost_up.objective_ = c_up_f;
    ranging.col_cost_dn.objective_ = c_dn_f;
    ranging.col_cost_up.in_var_ = c_up_e;
    ranging.col_cost_dn.in_var_ = c_dn_e;
    ranging.col_cost_up.ou_var_ = c_up_l;
    ranging.col_cost_dn.ou_var_ = c_dn_l;
  } else {
    // For maximization problems, flip data and negate the cost values
    ranging.col_cost_up.value_ = c_dn_c;
    ranging.col_cost_dn.value_ = c_up_c;
    ranging.col_cost_up.objective_ = c_dn_f;
    ranging.col_cost_dn.objective_ = c_up_f;
    ranging.col_cost_up.in_var_ = c_dn_e;
    ranging.col_cost_dn.in_var_ = c_up_e;
    ranging.col_cost_up.ou_var_ = c_dn_l;
    ranging.col_cost_dn.ou_var_ = c_up_l;
    std::transform(ranging.col_cost_up.value_.cbegin(),
                   ranging.col_cost_up.value_.cend(),
                   ranging.col_cost_up.value_.begin(), std::negate<double>());
    std::transform(ranging.col_cost_dn.value_.cbegin(),
                   ranging.col_cost_dn.value_.cend(),
                   ranging.col_cost_dn.value_.begin(), std::negate<double>());
  }

  ranging.col_bound_up.value_ = {b_up_b.begin(), b_up_b.begin() + numCol};
  ranging.col_bound_dn.value_ = {b_dn_b.begin(), b_dn_b.begin() + numCol};
  ranging.col_bound_up.objective_ = {b_up_f.begin(), b_up_f.begin() + numCol};
  ranging.col_bound_dn.objective_ = {b_dn_f.begin(), b_dn_f.begin() + numCol};
  ranging.col_bound_up.in_var_ = {b_up_e.begin(), b_up_e.begin() + numCol};
  ranging.col_bound_dn.in_var_ = {b_dn_e.begin(), b_dn_e.begin() + numCol};
  ranging.col_bound_up.ou_var_ = {b_up_l.begin(), b_up_l.begin() + numCol};
  ranging.col_bound_dn.ou_var_ = {b_dn_l.begin(), b_dn_l.begin() + numCol};

  // Flip all data and negate the row bound values
  ranging.row_bound_up.value_ = {b_dn_b.begin() + numCol,
                                 b_dn_b.begin() + numTotal};
  std::transform(ranging.row_bound_up.value_.cbegin(),
                 ranging.row_bound_up.value_.cend(),
                 ranging.row_bound_up.value_.begin(), std::negate<double>());

  ranging.row_bound_dn.value_ = {b_up_b.begin() + numCol,
                                 b_up_b.begin() + numTotal};
  std::transform(ranging.row_bound_dn.value_.cbegin(),
                 ranging.row_bound_dn.value_.cend(),
                 ranging.row_bound_dn.value_.begin(), std::negate<double>());

  ranging.row_bound_up.objective_ = {b_dn_f.begin() + numCol,
                                     b_dn_f.begin() + numTotal};
  ranging.row_bound_dn.objective_ = {b_up_f.begin() + numCol,
                                     b_up_f.begin() + numTotal};
  ranging.row_bound_up.in_var_ = {b_dn_e.begin() + numCol,
                                  b_dn_e.begin() + numTotal};
  ranging.row_bound_dn.in_var_ = {b_up_e.begin() + numCol,
                                  b_up_e.begin() + numTotal};
  ranging.row_bound_up.ou_var_ = {b_dn_l.begin() + numCol,
                                  b_dn_l.begin() + numTotal};
  ranging.row_bound_dn.ou_var_ = {b_up_l.begin() + numCol,
                                  b_up_l.begin() + numTotal};

  ranging.valid = true;
  if (solver_object.options_.log_dev_level)
    writeRangingFile(stdout, use_lp,
                     solver_object.highs_info_.objective_function_value,
                     solver_object.basis_, solver_object.solution_, ranging,
                     kSolutionStylePretty);
  return HighsStatus::kOk;
}

void writeRangingFile(FILE* file, const HighsLp& lp,
                      const double objective_function_value,
                      const HighsBasis& basis, const HighsSolution& solution,
                      const HighsRanging& ranging, const HighsInt style) {
  if (!ranging.valid) {
    fprintf(file, "None\n");
    return;
  }
  fprintf(file, "Valid\n");
  std::stringstream ss;
  assert(lp.col_names_.size() == static_cast<size_t>(lp.num_col_));
  assert(lp.row_names_.size() == static_cast<size_t>(lp.num_row_));
  const bool pretty = style == kSolutionStylePretty;
  const char* pretty_cost_format =
      "%6d   %4s  %-10.4g %-10.4g            %-10.4g            %-10.4g "
      "%-10.4g %-s\n";
  const char* raw_cost_format = "%-s %s %s %s %s\n";

  const char* pretty_bound_format =
      "%6d   %4s  %-10.4g %-10.4g %-10.4g %-10.4g %-10.4g %-10.4g %-10.4g "
      "%-s\n";
  const char* raw_bound_format = "%-s %s %s %s %s\n";

  std::array<char, 32> dn_objective;
  std::array<char, 32> up_objective;
  std::array<char, 32> dn_value;
  std::array<char, 32> up_value;

  auto objective = highsDoubleToString(objective_function_value,
                                       kRangingValueToStringTolerance);
  fprintf(file, "Objective %s\n", objective.data());
  if (pretty) {
    fprintf(file,
            "\n                                            Cost ranging\n"
            "Column Status  DownObj    Down                  Value             "
            "    Up         UpObj      Name\n");
  } else {
    fprintf(file, "\n# Cost ranging\n");
  }
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const std::string name = lp.col_names_[iCol];
    if (pretty) {
      fprintf(file, pretty_cost_format, (int)iCol,
              statusToString(basis.col_status[iCol], lp.col_lower_[iCol],
                             lp.col_upper_[iCol])
                  .c_str(),
              ranging.col_cost_dn.objective_[iCol],
              ranging.col_cost_dn.value_[iCol], lp.col_cost_[iCol],
              ranging.col_cost_up.value_[iCol],
              ranging.col_cost_up.objective_[iCol], name.c_str());
    } else {
      dn_objective = highsDoubleToString(ranging.col_cost_dn.objective_[iCol],
                                         kRangingValueToStringTolerance);
      up_objective = highsDoubleToString(ranging.col_cost_up.objective_[iCol],
                                         kRangingValueToStringTolerance);
      dn_value = highsDoubleToString(ranging.col_cost_dn.value_[iCol],
                                     kRangingValueToStringTolerance);
      up_value = highsDoubleToString(ranging.col_cost_up.value_[iCol],
                                     kRangingValueToStringTolerance);
      fprintf(file, raw_cost_format, name.c_str(), dn_objective.data(),
              dn_value.data(), up_value.data(), up_objective.data());
    }
  }
  if (pretty) {
    fprintf(file,
            "\n                                            Bound ranging\n"
            "Column Status  DownObj    Down       Lower      Value      Upper "
            "     Up         UpObj      Name\n");
  } else {
    fprintf(file, "\n# Bound ranging\n# Columns\n");
  }
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const std::string name = lp.col_names_[iCol];
    if (pretty) {
      fprintf(file, pretty_bound_format, (int)iCol,
              statusToString(basis.col_status[iCol], lp.col_lower_[iCol],
                             lp.col_upper_[iCol])
                  .c_str(),
              ranging.col_bound_dn.objective_[iCol],
              ranging.col_bound_dn.value_[iCol], lp.col_lower_[iCol],
              solution.col_value[iCol], lp.col_upper_[iCol],
              ranging.col_bound_up.value_[iCol],
              ranging.col_bound_up.objective_[iCol], name.c_str());
    } else {
      dn_objective = highsDoubleToString(ranging.col_bound_dn.objective_[iCol],
                                         kRangingValueToStringTolerance);
      up_objective = highsDoubleToString(ranging.col_bound_up.objective_[iCol],
                                         kRangingValueToStringTolerance);
      dn_value = highsDoubleToString(ranging.col_bound_dn.value_[iCol],
                                     kRangingValueToStringTolerance);
      up_value = highsDoubleToString(ranging.col_bound_up.value_[iCol],
                                     kRangingValueToStringTolerance);
      fprintf(file, raw_bound_format, name.c_str(), dn_objective.data(),
              dn_value.data(), up_value.data(), up_objective.data());
    }
  }

  if (pretty) {
    fprintf(file,
            "                                            Bound ranging\n"
            "   Row Status  DownObj    Down       Lower      Value      Upper "
            "     Up         UpObj      Name\n");
  } else {
    fprintf(file, "# Rows\n");
  }
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    const std::string name = lp.row_names_[iRow];
    if (pretty) {
      fprintf(file, pretty_bound_format, (int)iRow,
              statusToString(basis.row_status[iRow], lp.row_lower_[iRow],
                             lp.row_upper_[iRow])
                  .c_str(),
              ranging.row_bound_dn.objective_[iRow],
              ranging.row_bound_dn.value_[iRow], lp.row_lower_[iRow],
              solution.row_value[iRow], lp.row_upper_[iRow],
              ranging.row_bound_up.value_[iRow],
              ranging.row_bound_up.objective_[iRow], name.c_str());
    } else {
      dn_objective = highsDoubleToString(ranging.row_bound_dn.objective_[iRow],
                                         kRangingValueToStringTolerance);
      up_objective = highsDoubleToString(ranging.row_bound_up.objective_[iRow],
                                         kRangingValueToStringTolerance);
      dn_value = highsDoubleToString(ranging.row_bound_dn.value_[iRow],
                                     kRangingValueToStringTolerance);
      up_value = highsDoubleToString(ranging.row_bound_up.value_[iRow],
                                     kRangingValueToStringTolerance);
      fprintf(file, raw_bound_format, name.c_str(), dn_objective.data(),
              dn_value.data(), up_value.data(), up_objective.data());
    }
  }
}
