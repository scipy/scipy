/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsRanging.cpp
 * @brief Compute LP ranging data for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsRanging.h"

#include <algorithm>
#include <cassert>
#include <functional>  // for negate

#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"

double infProduct(double value) {
  // Multiplying value and HIGHS_CONST_INF
  if (value == 0) {
    return 0;
  } else {
    return value * HIGHS_CONST_INF;
  }
}

double possInfProduct(double poss_inf, double value) {
  // Multiplying something that could be infinite and value
  if (value == 0) {
    return 0;
  } else {
    return poss_inf * value;
  }
}

HighsStatus getHighsRanging(HighsRanging& ranging,
                            const HighsModelObject& highs_model_object) {
  if (highs_model_object.scaled_model_status_ != HighsModelStatus::OPTIMAL) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::ERROR,
                    "Cannot get ranging without an optimal solution");
    return HighsStatus::Error;
  }
  if (!highs_model_object.simplex_lp_status_.valid) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::ERROR,
                    "Cannot get ranging without a valid Simplex LP");
    return HighsStatus::Error;
  }
  // Aliases
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  const vector<double>& col_scale = highs_model_object.scale_.col_;
  const vector<double>& row_scale = highs_model_object.scale_.row_;
  const vector<double>& value_ = simplex_info.workValue_;
  const vector<double>& dual_ = simplex_info.workDual_;
  const vector<double>& cost_ = simplex_info.workCost_;
  const vector<double>& lower_ = simplex_info.workLower_;
  const vector<double>& upper_ = simplex_info.workUpper_;
  const vector<double>& Bvalue_ = simplex_info.baseValue_;
  const vector<double>& Blower_ = simplex_info.baseLower_;
  const vector<double>& Bupper_ = simplex_info.baseUpper_;
  const vector<int>& Nflag_ = simplex_basis.nonbasicFlag_;
  const vector<int>& Nmove_ = simplex_basis.nonbasicMove_;
  const vector<int>& Bindex_ = simplex_basis.basicIndex_;
  const HMatrix& matrix = highs_model_object.matrix_;
  const HFactor& factor = highs_model_object.factor_;

  // Local copies of scalars

  const int numRow = highs_model_object.simplex_lp_.numRow_;
  const int numCol = highs_model_object.simplex_lp_.numCol_;
  const int numTotal = numCol + numRow;
  const double H_TT = 1e-13;
  const double H_INF = HIGHS_CONST_INF;
  const double objective =
      highs_model_object.unscaled_solution_params_.objective_function_value;

  // Code written for minimization problems. Maximization problems are
  // solved by using negated costs in the simplex solver and
  // minimizing. Thus dual information in the simplex solver is
  // negated for maximization problems. The objective has the right
  // sign, though. Maximization problems are, thus, accommodated by
  // applying the sign multiplier to dual information.
  int sense = 1;
  if (highs_model_object.lp_.sense_ == ObjSense::MAXIMIZE) sense = -1;

  vector<int> iWork_(numTotal);
  vector<double> dWork_(numTotal);
  HVector column;
  column.setup(numRow);
  // As rgda/HModel.cpp HModel::sense() {

  vector<double> xi = Bvalue_;
  for (int i = 0; i < numRow; i++) {
    xi[i] = max(xi[i], Blower_[i]);
    xi[i] = min(xi[i], Bupper_[i]);
  }

  vector<double> dj = dual_;
  for (int j = 0; j < numTotal; j++) {
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
  for (int i = 0; i < numRow; i++) {
    dxi_inc[i] = Bupper_[i] - xi[i];
    dxi_dec[i] = Blower_[i] - xi[i];
  }

  vector<double> ddj_inc(numTotal);
  vector<double> ddj_dec(numTotal);
  for (int j = 0; j < numTotal; j++) {
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
  vector<int> ixj_inc(numTotal, -1);             // i-out
  vector<int> wxj_inc(numTotal, 0);              // which bound is limiting
  vector<int> jxj_inc(numTotal, -1);             // j = n(i), (with bound flip)

  vector<double> txj_dec(numTotal, -THETA_INF);
  vector<double> axj_dec(numTotal, 0);
  vector<int> ixj_dec(numTotal, -1);
  vector<int> wxj_dec(numTotal, 0);
  vector<int> jxj_dec(numTotal, -1);

  vector<double> tci_inc(numRow, +THETA_INF);  // theta
  vector<double> aci_inc(numRow, 0);           // alpha
  vector<int> jci_inc(numRow, -1);             // column index

  vector<double> tci_dec(numRow, -THETA_INF);
  vector<double> aci_dec(numRow, 0);
  vector<int> jci_dec(numRow, -1);

  // Major "theta" loop
  for (int j = 0; j < numTotal; j++) {
    // Skip basic column
    if (!Nflag_[j]) continue;

    // Form updated column
    column.clear();
    matrix.collect_aj(column, j, 1);
    factor.ftran(column, 0);
    int nWork = 0;
    for (int k = 0; k < column.count; k++) {
      int iRow = column.index[k];
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
    int myk_inc = -1;
    int myk_dec = -1;
    for (int k = 0; k < nWork; k++) {
      int i = iWork_[k];
      double alpha = dWork_[k];
      double theta_inc = (alpha < 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      double theta_dec = (alpha > 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      if (myt_inc > theta_inc) myt_inc = theta_inc, myk_inc = k;
      if (myt_dec < theta_dec) myt_dec = theta_dec, myk_dec = k;
    }

    if (myk_inc != -1) {
      int i = iWork_[myk_inc];
      double alpha = dWork_[myk_inc];
      ixj_inc[j] = i;
      axj_inc[j] = alpha;
      txj_inc[j] = (alpha < 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      wxj_inc[j] = (alpha < 0 ? +1 : -1);
    }

    if (myk_dec != -1) {
      int i = iWork_[myk_dec];
      double alpha = dWork_[myk_dec];
      ixj_dec[j] = i;
      axj_dec[j] = alpha;
      txj_dec[j] = (alpha > 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      wxj_dec[j] = (alpha > 0 ? +1 : -1);
    }

    // Accumulated dual ratio test
    double myd_inc = ddj_inc[j];
    double myd_dec = ddj_dec[j];
    for (int k = 0; k < nWork; k++) {
      int i = iWork_[k];
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
  for (int j = 0; j < numTotal; j++) {
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
  vector<int> c_up_e(numTotal), c_dn_e(numTotal);
  vector<int> c_up_l(numTotal), c_dn_l(numTotal);

  //
  // Ranging 2.1. non-basic cost ranging
  //
  //  const int check_col = 2951;
  for (int j = 0; j < numCol; j++) {
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

  for (int i = 0; i < numRow; i++) {
    if (Bindex_[i] < numCol) {
      // Primal variable and its sign
      int j = Bindex_[i], je;
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
  vector<int> b_up_e(numTotal), b_dn_e(numTotal);
  vector<int> b_up_l(numTotal), b_dn_l(numTotal);

  //
  // Ranging 3.1. non-basic bounds ranging
  //
  for (int j = 0; j < numTotal; j++) {
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
        int i = ixj_inc[j];
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
        assert(lower_[j] > -HIGHS_CONST_INF);
        b_up_f[j] = objective + sense * (upper_[j] - lower_[j]) * dualv;
        b_up_e[j] = j;
        b_up_l[j] = j;
      }

      // Decrease x_j
      if (ixj_dec[j] != -1) {
        int i = ixj_dec[j];
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
        assert(upper_[j] < HIGHS_CONST_INF);
        b_dn_f[j] = objective + sense * (lower_[j] - upper_[j]) * dualv;
        b_dn_e[j] = j;
        b_dn_l[j] = j;
      }
    }
  }

  //
  // Ranging 3.2. basic bounds ranging
  //
  for (int i = 0; i < numRow; i++) {
    for (int dir = -1; dir <= 1; dir += 2) {
      int j = Bindex_[i];
      double& newx = dir == -1 ? b_dn_b[j] : b_up_b[j];
      double& newf = dir == -1 ? b_dn_f[j] : b_up_f[j];
      int& j_enter = dir == -1 ? b_dn_e[j] : b_up_e[j];
      int& j_leave = dir == -1 ? b_dn_l[j] : b_up_l[j];

      int j_in = dir == -1 ? jci_inc[i] : jci_dec[i];
      double a_in = dir == -1 ? aci_inc[i] : aci_dec[i];
      if (j_in != -1) {
        int jmove = Nmove_[j_in];
        int i_out = jmove > 0 ? ixj_inc[j_in] : ixj_dec[j_in];
        int j_out = jmove > 0 ? jxj_inc[j_in] : jxj_dec[j_in];
        int w_out = jmove > 0 ? wxj_inc[j_in] : wxj_dec[j_in];
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
  // Ranging 4.1. Scale back
  //
  for (int j = 0; j < numCol; j++) {
    c_up_c[j] /= (c_up_c[j] == +H_INF) ? 1 : col_scale[j];
    c_dn_c[j] /= (c_dn_c[j] == -H_INF) ? 1 : col_scale[j];
    b_up_b[j] *= (b_up_b[j] == +H_INF) ? 1 : col_scale[j];
    b_dn_b[j] *= (b_dn_b[j] == +H_INF) ? 1 : col_scale[j];
  }
  for (int i = 0, j = numCol; i < numRow; i++, j++) {
    b_up_b[j] /= (b_up_b[j] == +H_INF) ? 1 : row_scale[i];
    b_dn_b[j] /= (b_dn_b[j] == +H_INF) ? 1 : row_scale[i];
  }

  //
  // Ranging 4.1.1 Trim small value to zero
  //
  for (int j = 0; j < numCol; j++) {
    if (fabs(c_up_c[j]) < H_TT) c_up_c[j] = 0;
    if (fabs(c_dn_c[j]) < H_TT) c_dn_c[j] = 0;
    if (fabs(b_up_b[j]) < H_TT) b_up_b[j] = 0;
    if (fabs(b_dn_b[j]) < H_TT) b_dn_b[j] = 0;
  }
  for (int i = 0, j = numCol; i < numRow; i++, j++) {
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

  return HighsStatus::OK;
}
