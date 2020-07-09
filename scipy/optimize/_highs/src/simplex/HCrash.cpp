/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HCrash.cpp
 * @brief Bixby and Maros-style crash for the HiGHS simplex solver
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HCrash.h"

#include <cassert>
#include <set>
#include <string>
#include <vector>

#include "lp_data/HConst.h"
#include "simplex/HSimplex.h"
#include "util/HighsSort.h"

using std::abs;
using std::cout;
using std::flush;
using std::max;
using std::string;
using std::vector;

void HCrash::crash(const int pass_crash_strategy) {
  crash_strategy = pass_crash_strategy;
  HighsLp& simplex_lp = workHMO.simplex_lp_;
  if (simplex_lp.numRow_ == 0) return;
  numRow = simplex_lp.numRow_;
  numCol = simplex_lp.numCol_;
  numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  if (crash_strategy == SIMPLEX_CRASH_STRATEGY_BASIC
#ifdef HiGHSDEV
      || crash_strategy == SIMPLEX_CRASH_STRATEGY_TEST_SING
#endif
  ) {
    // First and last variable types are the only types for basis and
    // test singularity crashes
    crsh_f_vr_ty = crsh_vr_ty_non_bc;
    crsh_l_vr_ty = crsh_vr_ty_bc;
    crsh_num_vr_ty = crsh_vr_ty_bc + 1;
    crsh_mn_pri_v = crsh_vr_ty_non_bc;
    crsh_mx_pri_v = crsh_vr_ty_bc;
    crsh_no_act_pri_v = crsh_mn_pri_v;
  } else {
    // First and last variable types are fixed and free for standard
    // crashes
    crsh_f_vr_ty = crsh_vr_ty_fx;
    crsh_l_vr_ty = crsh_vr_ty_fr;
    crsh_num_vr_ty = crsh_vr_ty_fr + 1;
    crsh_mn_pri_v = crsh_vr_ty_fx;
    crsh_mx_pri_v = crsh_vr_ty_fr;
    crsh_no_act_pri_v = crsh_mn_pri_v;
  }

  if (crash_strategy == SIMPLEX_CRASH_STRATEGY_BIXBY ||
      crash_strategy == SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS) {
    // Use the Bixby crash
    bixby();
  }
#ifdef HiGHSDEV
  else if (crash_strategy == SIMPLEX_CRASH_STRATEGY_TEST_SING) {
    // Use the test singularity crash
    tsSing();
  }
#endif
  else {
    // Use the LTSSF crash
    ltssf();
  }
}

void HCrash::bixby() {
  HighsLp& simplex_lp = workHMO.simplex_lp_;

  const int* Astart = &simplex_lp.Astart_[0];
  const int* Aindex = &simplex_lp.Aindex_[0];
  const double* Avalue = &simplex_lp.Avalue_[0];

  bixby_no_nz_c_co =
      crash_strategy == SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS;
  bixby_no_nz_c_co = false;

  bool perform_crash = bixby_iz_da();
  if (!perform_crash) return;

  // bixby_rp_mrt(workHMO);

  // These multipliers are in Step 2(a) and Step 2(b) of the paper: default
  // values 0.99 and 0.01
  bixby_mu_a = 0.99;
  bixby_mu_b = 0.01;

#ifdef HiGHSDEV
  printf("\nBixby Crash");
  if (bixby_no_nz_c_co) {
    printf(": No basic columns with nonzero costs\n");
  } else {
    printf(": Any basic columns regardless of cost\n");
  }
#endif
  for (int ps_n = 0; ps_n < numCol; ps_n++) {
    //  In each pass:
    //  Consider column c_n
    int c_n = bixby_mrt_ix[ps_n];
    double c_mx_abs_v = crsh_mtx_c_mx_abs_v[c_n];
    //  Find the max |entry| over all rows with a count of zero
    int r_o_mx_aa = -1;
    double aa = 0;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      if (bixby_r_k[r_n] == 0) {
        double lc_aa = fabs(Avalue[el_n]);
        if (lc_aa > aa) {
          aa = lc_aa;
          r_o_mx_aa = r_n;
        }
      }
    }
#ifdef HiGHSDEV
    if (reportBixbyPass)
      printf("Pass %3d: c_n = %3d; MxEn = %10g; aa = %10g; aa/MxEn = %10g",
             ps_n, c_n, c_mx_abs_v, aa, aa / c_mx_abs_v);
#endif
    // Scale aa by the max |entry| in the column since CPLEX assumes
    // the matrix is scaled so the max entry in each column is 1
    aa /= c_mx_abs_v;
    bool nx_ps = false;
    if (aa >= bixby_mu_a) {
      assert(r_o_mx_aa >= 0);
      // Column pv_c_n becomes basic in row pv_r_n
      int pv_c_n = c_n;
      int pv_r_n = r_o_mx_aa;
      // printf(" ** Type a: c_n = %d; pv_c_n = %d ** c_n, pv_c_n);
      bixby_pv_in_r[pv_r_n] = 1;
      bixby_vr_in_r[pv_r_n] = pv_c_n;
      bixby_pseudo_pv_v[pv_r_n] = aa;
      for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
        // int r_n = Aindex[el_n];
        // printf("\n Row %3d: value %g", r_n, Avalue[el_n]);
        bixby_r_k[Aindex[el_n]] += 1;
      }
      bixby_n_cdd_r -= 1;
#ifdef HiGHSDEV
      if (reportBixbyPass)
        printf(": pv_r = %3d; n_cdd_r = %3d\n", pv_r_n, bixby_n_cdd_r);
#endif
      nx_ps = true;
    } else {
      // Find out if there is some row l for which |entry| > bixby_mu_b * v_l
#ifdef HiGHSDEV
      double rp_v;
#endif
      for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
        int r_n = Aindex[el_n];
        // If this value in the column would give an unacceptable
        // multiplier then continue to next pass
        nx_ps = fabs(Avalue[el_n]) >
                bixby_mu_b * bixby_pseudo_pv_v[r_n] * c_mx_abs_v;
        if (nx_ps) {
#ifdef HiGHSDEV
          rp_v = fabs(Avalue[el_n]) / (bixby_pseudo_pv_v[r_n] * c_mx_abs_v);
#endif
          break;
        }
      }
#ifdef HiGHSDEV
      if (nx_ps && reportBixbyPass)
        printf(": Unacceptable multiplier of %g > %g\n", rp_v, bixby_mu_b);
#endif
    }
    // Some value in the column would give an unacceptable multiplier
    // so continue to next pass
    if (nx_ps) continue;
    // Find out whether there is an entry in a row with no pivot
    aa = 0;
    r_o_mx_aa = no_ix;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      if (bixby_pv_in_r[r_n] == 0) {
        double lc_aa = fabs(Avalue[el_n]);
        if (lc_aa > aa) {
          aa = lc_aa;
          r_o_mx_aa = r_n;
        }
      }
    }
    // If there is no entry in a row with no pivot then continue to
    // next pass
    if (r_o_mx_aa == no_ix) {
#ifdef HiGHSDEV
      if (reportBixbyPass) printf(": No entry in a row with no pivot\n");
#endif
      continue;
    }
    // Scale aa by the max |entry| in the column since CPLEX assumes
    // the matrix is scaled so the max entry in each column is 1
    aa /= c_mx_abs_v;
    // Column pv_c_n becomes basic in row pv_r_n
    int pv_c_n = c_n;
    int pv_r_n = r_o_mx_aa;
    bixby_pv_in_r[pv_r_n] = 1;
    bixby_vr_in_r[pv_r_n] = pv_c_n;
    bixby_pseudo_pv_v[pv_r_n] = aa;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      bixby_r_k[Aindex[el_n]] += 1;
    }
    bixby_n_cdd_r -= 1;
#ifdef HiGHSDEV
    if (reportBixbyPass)
      printf(": pv_r = %3d; n_cdd_r = %3d\n", pv_r_n, bixby_n_cdd_r);
#endif
    if (bixby_n_cdd_r == 0) break;
  }
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (bixby_vr_in_r[r_n] == no_ix) continue;
    if (bixby_vr_in_r[r_n] == numCol + r_n) continue;
    int cz_c_n = bixby_vr_in_r[r_n];
    int columnIn = cz_c_n;
    int columnOut = numCol + r_n;
    workHMO.simplex_basis_.nonbasicFlag_[columnIn] = NONBASIC_FLAG_FALSE;
    workHMO.simplex_basis_.nonbasicFlag_[columnOut] = NONBASIC_FLAG_TRUE;
#ifdef HiGHSDEV
    int cz_r_n = r_n;
    int vr_ty = crsh_r_ty[cz_r_n];
    crsh_vr_ty_rm_n_r[vr_ty] += 1;
    vr_ty = crsh_c_ty[cz_c_n];
    crsh_vr_ty_add_n_c[vr_ty] += 1;
#endif
  }
#ifdef HiGHSDEV
  // Analyse the row and column status after Crash
  // basicIndex is only required for this analysis, so set it here.
  workHMO.simplex_basis_.basicIndex_.resize(numRow);
  initialise_basic_index(workHMO);
  crsh_an_r_c_st_af();
#endif
}

bool HCrash::bixby_iz_da() {
  HighsLp& simplex_lp = workHMO.simplex_lp_;
  const int* Astart = &simplex_lp.Astart_[0];
  const double* Avalue = &simplex_lp.Avalue_[0];
  const int objSense = (int)simplex_lp.sense_;
  const double* colCost = &simplex_lp.colCost_[0];
  const double* colLower = &simplex_lp.colLower_[0];
  const double* colUpper = &simplex_lp.colUpper_[0];

  // const double *primalColLowerImplied = simplex_lp.primalColLowerImplied_;
  // const double *primalColUpperImplied = simplex_lp.primalColUpperImplied_;
  //
  // const double *dualColLowerImplied = simplex_lp.dualColLowerImplied_;
  // const double *dualColUpperImplied = simplex_lp.dualColUpperImplied_;

  // Allocate the arrays required for crash
  crsh_mtx_c_mx_abs_v.resize(numCol);

  bixby_mrt_v.resize(numCol + 1);
  bixby_pseudo_pv_v.resize(numRow);
  bixby_mrt_ix.resize(numCol + 1);
  bixby_pv_in_r.resize(numRow);
  bixby_vr_in_r.resize(numRow);
  bixby_r_k.resize(numRow);
  // bixby_ze_r_k.resize(numRow);

#ifdef HiGHSDEV
  crsh_an_c_co();
#endif
  crsh_iz_vr_ty();

#ifdef HiGHSDEV
  crsh_rp_r_c_st(0);
#endif

  // Initialise the arrays required for the Bixby crash
  //
  // bixby_pseudo_pv_v: "v" in the paper: this is a pseudo pivot value
  // for each row
  //
  // bixby_pv_in_r: "I" in the paper: this is a 0/1 flag to indicate
  // whether there is a basic variable in each row
  //
  // bixby_vr_in_r: "B" is the paper: this is the basic variable in a
  // particular row
  //
  // bixby_r_k: "r" in the paper: this is the number of entries in
  // each row of the basis matrix
  bixby_n_cdd_r = numRow;
  for (int r_n = 0; r_n < numRow; r_n++) {
    bixby_pseudo_pv_v[r_n] = HIGHS_CONST_INF;
    if (crsh_r_ty[r_n] == crsh_vr_ty_fx) {
      bixby_pv_in_r[r_n] = 0;
      bixby_vr_in_r[r_n] = no_ix;
      bixby_r_k[r_n] = 0;
      // bixby_ze_r_k[r_n] = 1;
    } else {
      bixby_pv_in_r[r_n] = 1;
      bixby_vr_in_r[r_n] = numCol + r_n;
      bixby_r_k[r_n] = 1;
      // bixby_ze_r_k[r_n] = 0;
      bixby_n_cdd_r -= 1;
    }
  }
  if (bixby_n_cdd_r == 0) return false;
  double mx_co_v = -HIGHS_CONST_INF;
  for (int c_n = 0; c_n < numCol; c_n++) {
    // Find largest |entry| in each column
    crsh_mtx_c_mx_abs_v[c_n] = 0.0;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      crsh_mtx_c_mx_abs_v[c_n] =
          max(fabs(Avalue[el_n]), crsh_mtx_c_mx_abs_v[c_n]);
    }
    double sense_col_cost = objSense * colCost[c_n];
    mx_co_v = max(fabs(sense_col_cost), mx_co_v);
  }
  double co_v_mu = 1;
  if (mx_co_v > 0) co_v_mu = 1e3 * mx_co_v;
  // ... and then updated with c_j/c_max, being colCost[c_n]/co_v_mu
  // So, first compute the cost coefficient of maximum absolute value
  int os;
  int n_en;
  os = 0;
  // Free columns - impossible after presolve
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_fr) continue;
    double sense_col_cost = objSense * colCost[c_n];
    if (bixby_no_nz_c_co && sense_col_cost != 0.0) continue;
    n_en += 1;
    bixby_mrt_v[os + n_en] = sense_col_cost / co_v_mu;
    bixby_mrt_ix[os + n_en] = c_n;
  }
  if (n_en > 0) {
    maxheapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os + 1; en <= os + n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en + 1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en + 1];
    }
    os += 1 + n_en;
  }
#ifdef HiGHSDEV
  printf("%5d Free    cols (%1d)\n", n_en, crsh_vr_ty_fr);
#endif

  // 1-sided columns
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_1_sd) continue;
    double sense_col_cost = objSense * colCost[c_n];
    if (bixby_no_nz_c_co && sense_col_cost != 0.0) continue;
    n_en += 1;
    if (colUpper[c_n] >= HIGHS_CONST_INF) {
      bixby_mrt_v[os + n_en] = colLower[c_n] + sense_col_cost / co_v_mu;
    } else {
      bixby_mrt_v[os + n_en] = -colUpper[c_n] + sense_col_cost / co_v_mu;
    }
    bixby_mrt_ix[os + n_en] = c_n;
  }
  if (n_en > 0) {
    maxheapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os + 1; en <= os + n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en + 1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en + 1];
    }
    os += 1 + n_en;
  }
#ifdef HiGHSDEV
  printf("%5d 1-sided cols (%1d)\n", n_en, crsh_vr_ty_1_sd);
#endif

  // 2-sided columns
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_2_sd) continue;
    double sense_col_cost = objSense * colCost[c_n];
    if (bixby_no_nz_c_co && sense_col_cost != 0.0) continue;
    n_en += 1;
    bixby_mrt_v[os + n_en] =
        colLower[c_n] - colUpper[c_n] + sense_col_cost / co_v_mu;
    bixby_mrt_ix[os + n_en] = c_n;
  }
  if (n_en > 0) {
    maxheapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os + 1; en <= os + n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en + 1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en + 1];
    }
    os += 1 + n_en;
  }
#ifdef HiGHSDEV
  printf("%5d 2-sided cols (%1d)\n", n_en, crsh_vr_ty_2_sd);
#endif

  // Fixed columns - impossible after presolve
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_fx) continue;
    double sense_col_cost = objSense * colCost[c_n];
    if (bixby_no_nz_c_co && sense_col_cost != 0.0) continue;
    n_en += 1;
    bixby_mrt_v[os + n_en] = HIGHS_CONST_INF;
    bixby_mrt_ix[os + n_en] = c_n;
  }
  if (n_en > 0) {
    maxheapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os + 1; en <= os + n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en + 1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en + 1];
    }
    os += 1 + n_en;
  }
#ifdef HiGHSDEV
  printf("%5d Fixed   cols (%1d)\n", n_en, crsh_vr_ty_fx);
#endif
  return true;
}

void HCrash::bixby_rp_mrt() {
  HighsLp& simplex_lp = workHMO.simplex_lp_;
  const int objSense = (int)simplex_lp.sense_;
  const double* colCost = &simplex_lp.colCost_[0];
  const double* colLower = &simplex_lp.colLower_[0];
  const double* colUpper = &simplex_lp.colUpper_[0];
  double mx_co_v = -HIGHS_CONST_INF;
  for (int c_n = 0; c_n < numCol; c_n++) {
    double sense_col_cost = objSense * colCost[c_n];
    mx_co_v = max(fabs(sense_col_cost), mx_co_v);
  }
  double co_v_mu = 1;
  if (mx_co_v > 0) co_v_mu = 1e3 * mx_co_v;
  double prev_mrt_v0 = -HIGHS_CONST_INF;
  double prev_mrt_v = -HIGHS_CONST_INF;
  bool rp_c;
  bool rp_al_c = false;
  int n_mrt_v = 0;
  if (mx_co_v > 0) co_v_mu = 1e3 * mx_co_v;
  printf("\nAnalysis of sorted Bixby merits\n");
  for (int ps_n = 0; ps_n < numCol; ps_n++) {
    double mrt_v = bixby_mrt_v[ps_n];
    int c_n = bixby_mrt_ix[ps_n];
    double sense_col_cost = objSense * colCost[c_n];
    double mrt_v0 = mrt_v - sense_col_cost / co_v_mu;
    double c_lb = colLower[c_n];
    double c_ub = colUpper[c_n];
    if ((ps_n == 0) || (ps_n == numCol - 1))
      rp_c = true;
    else if ((crsh_c_ty[c_n] != crsh_c_ty[bixby_mrt_ix[ps_n - 1]]) ||
             (crsh_c_ty[c_n] != crsh_c_ty[bixby_mrt_ix[ps_n + 1]])) {
      rp_c = true;
      prev_mrt_v = -HIGHS_CONST_INF;
      prev_mrt_v0 = -HIGHS_CONST_INF;
    } else if (rp_al_c)
      rp_c = true;
    else
      rp_c = mrt_v0 > prev_mrt_v0;
    prev_mrt_v0 = mrt_v0;
    if (mrt_v > prev_mrt_v) {
      n_mrt_v += 1;
      prev_mrt_v = mrt_v;
    }
    if (rp_c)
      printf(
          "%5d: Col %5d, Type = %1d; MrtV = %10.4g; MrtV0 = %10.4g; "
          "[%10.4g,%10.4g]\n",
          ps_n, c_n, crsh_c_ty[c_n], mrt_v, mrt_v0, c_lb, c_ub);
  }
  printf("\n%6d different Bixby merits\n", n_mrt_v);
}

void HCrash::ltssf() {
  HighsLp& simplex_lp = workHMO.simplex_lp_;
  if (crash_strategy == SIMPLEX_CRASH_STRATEGY_LTSSF_K) {
    crsh_fn_cf_pri_v = 1;
    crsh_fn_cf_k = 10;
    alw_al_bs_cg = false;
    no_ck_pv = false;
  } else if (crash_strategy == SIMPLEX_CRASH_STRATEGY_LTSF_K) {
    crsh_fn_cf_pri_v = 1;
    crsh_fn_cf_k = 10;
    alw_al_bs_cg = false;
    no_ck_pv = true;
  } else if (crash_strategy == SIMPLEX_CRASH_STRATEGY_LTSF) {
    crsh_fn_cf_pri_v = 1;
    crsh_fn_cf_k = 10;
    alw_al_bs_cg = true;
    no_ck_pv = true;
  } else if (crash_strategy == SIMPLEX_CRASH_STRATEGY_LTSSF_PRI) {
    crsh_fn_cf_pri_v = 10;
    crsh_fn_cf_k = 1;
    alw_al_bs_cg = false;
    no_ck_pv = false;
  } else if (crash_strategy == SIMPLEX_CRASH_STRATEGY_LTSF_PRI) {
    crsh_fn_cf_pri_v = 10;
    crsh_fn_cf_k = 1;
    alw_al_bs_cg = false;
    no_ck_pv = false;
  } else if (crash_strategy == SIMPLEX_CRASH_STRATEGY_BASIC) {
    crsh_fn_cf_pri_v = 10;
    crsh_fn_cf_k = 1;
    alw_al_bs_cg = false;
    no_ck_pv = false;
  } else {
    //  Dev version
    crsh_fn_cf_pri_v = 1;
    crsh_fn_cf_k = 10;
    alw_al_bs_cg = false;
    no_ck_pv = false;
  }

  mn_co_tie_bk = false;
  numRow = simplex_lp.numRow_;
  numCol = simplex_lp.numCol_;
  numTot = simplex_lp.numCol_ + simplex_lp.numRow_;

  // Initialise the LTSSF data structures
  ltssf_iz_da();
#ifdef HiGHSDEV
  printf("\nLTSSF Crash\n");
  printf(" crsh_fn_cf_pri_v = %d\n", crsh_fn_cf_pri_v);
  printf(" crsh_fn_cf_k = %d\n", crsh_fn_cf_k);
  printf(" Objective f = %2d*pri_v - %2d*k\n", crsh_fn_cf_pri_v, crsh_fn_cf_k);
  string lg;
  // if (mn_co_tie_bk) lg = "T"; else lg = "F";
  // printf(" Break priority function ties on smallest cost = %s\n", lg);
  if (alw_al_bs_cg)
    printf(" Allow any basis change regardless of priority = T\n");
  else
    printf(" Allow any basis change regardless of priority = F\n");
  if (no_ck_pv)
    printf(" Don't do numerical check on pivot = T\n");
  else
    printf(" Don't do numerical check on pivot = F\n");

  printf("Max row priority is %d\n", mx_r_pri);
  printf("Max col priority is %d\n", mx_c_pri);
#endif
  if ((!alw_al_bs_cg) && (mx_r_pri + mx_c_pri <= crsh_mx_pri_v)) {
#ifdef HiGHSDEV
    printf(
        "Max row priority of %d + Max col priority of %d = %d <= %d: no value "
        "in performing LTSSF crash\n",
        mx_r_pri, mx_c_pri, mx_r_pri + mx_c_pri, crsh_mx_pri_v);
#endif
    // Save the solved results
    return;
  }
#ifdef HiGHSDEV
  if (ltssf_ck_fq > 0) {
    printf("\nCHECKING LTSSF DATA NOW AND EVERY %d PASS(ES)!!\n\n",
           ltssf_ck_fq);
    ltssf_ck_da();
  }
  if (reportCrashData) ltssf_rp_pri_k_da();
#endif

#ifdef HiGHSDEV
  crsh_rp_r_c_st(0);
#endif
  ltssf_iterate();

#ifdef HiGHSDEV
  printf(" %d/%d basis changes from %d passes\n", n_crsh_bs_cg, numRow,
         n_crsh_ps);
  printf(
      " Absolute tolerance (%6.4f): Rejected %7d pivots: min absolute pivot "
      "value = %6.4e\n",
      tl_crsh_abs_pv_v, n_abs_pv_no_ok, mn_abs_pv_v);
  printf(
      " Relative tolerance (%6.4f): Rejected %7d pivots: min relative pivot "
      "value = %6.4e\n",
      tl_crsh_rlv_pv_v, n_rlv_pv_no_ok, mn_rlv_pv_v);
  // Analyse the row and column status after Crash
  // basicIndex is only required for this analysis, so set it here.
  workHMO.simplex_basis_.basicIndex_.resize(numRow);
  initialise_basic_index(workHMO);
  crsh_an_r_c_st_af();
#endif
}

void HCrash::ltssf_iz_mode() {
  crsh_fn_cf_pri_v = 1;
  crsh_fn_cf_k = 10;
  alw_al_bs_cg = false;
  no_ck_pv = false;
}

void HCrash::ltssf_iterate() {
  // LTSSF Main loop
  n_crsh_ps = 0;
  n_crsh_bs_cg = 0;
  bool ltssf_stop = false;
  for (;;) {
    ltssf_cz_r();
    if (cz_r_n == no_ix) break;
    cz_r_pri_v = crsh_r_ty_pri_v[crsh_r_ty[cz_r_n]];
    ltssf_cz_c();
    bool bs_cg = cz_c_n != no_ix;
    if (bs_cg) {
#ifdef HiGHSDEV
      if (reportCrashData)
        printf("Pass %2d; cz_r = %2d; cz_c = %2d\n", n_crsh_ps, cz_r_n, cz_c_n);
#endif
      // A basis change has occurred
      n_crsh_bs_cg += 1;
      double abs_pv_v = fabs(pv_v);
      double rlv_pv_v = abs_pv_v / crsh_mtx_c_mx_abs_v[cz_c_n];
      mn_abs_pv_v = min(abs_pv_v, mn_abs_pv_v);
      mn_rlv_pv_v = min(rlv_pv_v, mn_rlv_pv_v);
      int columnIn = cz_c_n;
      int columnOut = numCol + cz_r_n;
      workHMO.simplex_basis_.nonbasicFlag_[columnIn] = NONBASIC_FLAG_FALSE;
      workHMO.simplex_basis_.nonbasicFlag_[columnOut] = NONBASIC_FLAG_TRUE;
      // Update the count of this type of removal and addition
#ifdef HiGHSDEV
      int vr_ty = crsh_r_ty[cz_r_n];
      crsh_vr_ty_rm_n_r[vr_ty] += 1;
      vr_ty = crsh_c_ty[cz_c_n];
      crsh_vr_ty_add_n_c[vr_ty] += 1;
#endif
    } else {
#ifdef HiGHSDEV
      if (reportCrashData)
        printf("Pass %2d; cz_r = %2d: No basis change\n", n_crsh_ps, cz_r_n);
#endif
    }
#ifdef HiGHSDEV
    if (reportCrashData) ltssf_rp_pri_k_da();
#endif
    ltssf_u_da();
    // Check LTSSF data every ltssf_ck_fq passes (if ltssf_ck_fq>0)
#ifdef HiGHSDEV
    if ((ltssf_ck_fq > 0) && (n_crsh_ps % ltssf_ck_fq == 0)) ltssf_ck_da();
    if (reportCrashData) ltssf_rp_pri_k_da();
#endif
    // Determine whether the are still rows worth removing
    mx_r_pri = crsh_mn_pri_v - 1;
    for (int pri_v = crsh_mx_pri_v; pri_v > crsh_mn_pri_v; pri_v--) {
      if (crsh_r_pri_mn_r_k[pri_v] < numCol + 1) {
        mx_r_pri = pri_v;
        break;
      }
    }
    if ((!alw_al_bs_cg) && (mx_r_pri + mx_c_pri <= crsh_mx_pri_v)) {
#ifdef HiGHSDEV
      printf(
          "Max active row priority of %d + Max original col priority of %d = "
          "%d <= %d: no value in performing further LTSSF crash\n",
          mx_r_pri, mx_c_pri, mx_r_pri + mx_c_pri, crsh_mx_pri_v);
#endif
      ltssf_stop = true;
    }
    n_crsh_ps += 1;
    if (ltssf_stop) break;
  }
}

void HCrash::ltssf_u_da() {
  if ((cz_r_n != no_ix) && (cz_c_n != no_ix)) {
    ltssf_u_da_af_bs_cg();
  } else {
    ltssf_u_da_af_no_bs_cg();
  }
  // If there are no more rows with the current maximum row priority
  // then get the new maximum row priority value TODO Surely this is
  // not necessary with 2-d headers
  if ((crsh_r_pri_mn_r_k[cz_r_pri_v] > numCol) && (cz_r_pri_v == mx_r_pri_v)) {
    mx_r_pri_v = -HIGHS_CONST_I_INF;
    for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++)
      if (crsh_r_pri_mn_r_k[pri_v] <= numCol) mx_r_pri_v = pri_v;
  }
}

void HCrash::ltssf_u_da_af_bs_cg() {
  HighsLp& simplex_lp = workHMO.simplex_lp_;
  const int* Astart = &simplex_lp.Astart_[0];
  const int* Aindex = &simplex_lp.Aindex_[0];
  // ltssf_rp_r_k();
  for (int r_el_n = CrshARstart[cz_r_n]; r_el_n < CrshARstart[cz_r_n + 1];
       r_el_n++) {
    int c_n = CrshARindex[r_el_n];
    if (crsh_act_c[c_n] == crsh_vr_st_no_act) continue;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      if (crsh_act_r[r_n] == crsh_vr_st_no_act) continue;
      // Remove the row from the linked list with this number of active entries
      int prev_r_n;
      int r_k = crsh_r_k[r_n];
      int pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
#ifdef HiGHSDEV
      if (reportCrashData) {
        ltssf_rp_pri_k_da();
        printf("1: Remove row %d of pri %d from linked list with %d entries\n",
               r_n, pri_v, r_k);
      }
#endif
      int hdr_ix = pri_v * (numCol + 1) + r_k;
      // Remove the row from the linked list with this number of active entries
      int nx_r_n = crsh_r_pri_k_lkf[r_n];
      if (r_n == crsh_r_pri_k_hdr[hdr_ix]) {
        prev_r_n = no_lk;
        crsh_r_pri_k_hdr[hdr_ix] = nx_r_n;
      } else {
        prev_r_n = crsh_r_pri_k_lkb[r_n];
        crsh_r_pri_k_lkf[prev_r_n] = nx_r_n;
      }
      if (nx_r_n != no_lk) crsh_r_pri_k_lkb[nx_r_n] = prev_r_n;
      if ((crsh_r_pri_k_hdr[hdr_ix] == no_lk) &&
          (crsh_r_pri_mn_r_k[pri_v] == r_k)) {
        // This was the only row of minimum row count so look for the next row
        // count with non-null header
        //
        // Set crsh_r_pri_mn_r_k to numCol+1 in case r_k=numCol so priority is
        // cleared
        crsh_r_pri_mn_r_k[pri_v] = numCol + 1;
        for (int qy_k = r_k + 1; qy_k < numCol + 1; qy_k++) {
          int hdr_ix = pri_v * (numCol + 1) + qy_k;
          if (crsh_r_pri_k_hdr[hdr_ix] != no_lk) {
            crsh_r_pri_mn_r_k[pri_v] = qy_k;
            break;
          }
        }
      }
      // Reduce the number of active entries in this row by one and...
      r_k -= 1;
      crsh_r_k[r_n] = r_k;
      if (r_k > 0) {
        // ... either add the row as the header of the list with one
        // fewer number of active entries...
#ifdef HiGHSDEV
        if (reportCrashData) {
          ltssf_rp_pri_k_da();
          printf("Add row %d of pri %d to linked list with %d entries\n", r_n,
                 pri_v, r_k);
        }
#endif
        int hdr_ix = pri_v * (numCol + 1) + r_k;
        nx_r_n = crsh_r_pri_k_hdr[hdr_ix];
        crsh_r_pri_k_hdr[hdr_ix] = r_n;
        crsh_r_pri_k_lkf[r_n] = nx_r_n;
        if (nx_r_n != no_lk) crsh_r_pri_k_lkb[nx_r_n] = r_n;
        if (crsh_r_pri_mn_r_k[pri_v] > r_k) {
          // There is now a row of smaller count for this priority
          crsh_r_pri_mn_r_k[pri_v] = r_k;
        }
      } else {
#ifdef HiGHSDEV
        if (reportCrashData) {
          ltssf_rp_pri_k_da();
          printf(
              "2: Remove row %d of pri %d and count %d from active submatrix\n",
              r_n, pri_v, r_k);
        }
#endif
        // ...or, if the count is zero, the row leaves the active submatrix...
        crsh_act_r[r_n] = crsh_vr_st_no_act;
        // ... and has already left the priority value and count data structure
      }
    }
    // The column leaves the active submatrix
    crsh_act_c[c_n] = crsh_vr_st_no_act;
    // Unless this is the chosen column, the structural remains nonbasic.
  }
}

void HCrash::ltssf_u_da_af_no_bs_cg() {
  // The basis has not changed. The row becomes inactive and the
  // number of active entries in each column in which it has entries
  // is reduced by one.
  for (int r_el_n = CrshARstart[cz_r_n]; r_el_n < CrshARstart[cz_r_n + 1];
       r_el_n++) {
    int c_n = CrshARindex[r_el_n];
    if (crsh_act_c[c_n] == crsh_vr_st_no_act) continue;
    crsh_c_k[c_n] -= 1;
    // If the number of active entries in the column is zeroed then it becomes
    // inactive.
    if (crsh_c_k[c_n] == 0) crsh_act_c[c_n] = crsh_vr_st_no_act;
  }
  int r_n = cz_r_n;
  // Remove the row from the linked list with this number of active entries
  // Remove the row from the linked list with this priority and number of active
  // entries
  crsh_act_r[r_n] = crsh_vr_st_no_act;
  int pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
  int r_k = crsh_r_k[r_n];
  int hdr_ix = pri_v * (numCol + 1) + r_k;
  int prev_r_n;
  int nx_r_n = crsh_r_pri_k_lkf[r_n];
  if (r_n == crsh_r_pri_k_hdr[hdr_ix]) {
    prev_r_n = no_lk;
    crsh_r_pri_k_hdr[hdr_ix] = nx_r_n;
  } else {
    prev_r_n = crsh_r_pri_k_lkb[r_n];
    crsh_r_pri_k_lkf[prev_r_n] = nx_r_n;
  }
  if (nx_r_n != no_lk) crsh_r_pri_k_lkb[nx_r_n] = prev_r_n;
  if ((crsh_r_pri_k_hdr[hdr_ix] == no_lk) &&
      (crsh_r_pri_mn_r_k[pri_v] == r_k)) {
    // No more rows of this count for this priority value: look for next
    // highest.
    //
    // Set crsh_r_pri_mn_r_k to numCol + 1 in case r_k=numCol so
    // priority is cleared
    crsh_r_pri_mn_r_k[pri_v] = numCol + 1;
    for (int qy_k = r_k + 1; qy_k < numCol + 1; qy_k++) {
      int hdr_ix = pri_v * (numCol + 1) + qy_k;
      if (crsh_r_pri_k_hdr[hdr_ix] != no_lk) {
        crsh_r_pri_mn_r_k[pri_v] = qy_k;
        break;
      }
    }
  }
}

void HCrash::ltssf_iz_da() {
  HighsLp& simplex_lp = workHMO.simplex_lp_;
  SimplexBasis& simplex_basis = workHMO.simplex_basis_;
  // bool ImpliedDualLTSSF = false;
  // ImpliedDualLTSSF = true;
  const int* Astart = &simplex_lp.Astart_[0];
  const int* Aindex = &simplex_lp.Aindex_[0];
  const double* Avalue = &simplex_lp.Avalue_[0];
  int numEl = Astart[numCol];
  // const double *primalColLowerImplied = simplex_lp.primalColLowerImplied_;
  // const double *primalColUpperImplied = simplex_lp.primalColUpperImplied_;
  // const double *primalRowLowerImplied = simplex_lp.primalRowLowerImplied_;
  // const double *primalRowUpperImplied = simplex_lp.primalRowUpperImplied_;
  //
  // const double *dualColLowerImplied = simplex_lp.dualColLowerImplied_;
  // const double *dualColUpperImplied = simplex_lp.dualColUpperImplied_;
  // const double *dualRowLowerImplied = simplex_lp.dualRowLowerImplied_;
  // const double *dualRowUpperImplied = simplex_lp.dualRowUpperImplied_;

  // Allocate the crash variable type arrays
  crsh_r_ty_pri_v.resize(crsh_num_vr_ty);
  crsh_c_ty_pri_v.resize(crsh_num_vr_ty);
  if (crash_strategy == SIMPLEX_CRASH_STRATEGY_BASIC) {
    // Basis-preserving crash:
    crsh_r_ty_pri_v[crsh_vr_ty_non_bc] = 1;
    crsh_r_ty_pri_v[crsh_vr_ty_bc] = 0;
    crsh_c_ty_pri_v[crsh_vr_ty_non_bc] = 0;
    crsh_c_ty_pri_v[crsh_vr_ty_bc] = 1;

  } else {
    // Standard crash:
    crsh_r_ty_pri_v[crsh_vr_ty_fx] = 3;
    crsh_r_ty_pri_v[crsh_vr_ty_2_sd] = 2;
    crsh_r_ty_pri_v[crsh_vr_ty_1_sd] = 1;
    crsh_r_ty_pri_v[crsh_vr_ty_fr] = 0;
    crsh_c_ty_pri_v[crsh_vr_ty_fx] = 0;
    crsh_c_ty_pri_v[crsh_vr_ty_2_sd] = 1;
    crsh_c_ty_pri_v[crsh_vr_ty_1_sd] = 2;
    crsh_c_ty_pri_v[crsh_vr_ty_fr] = 3;
  }

  // Allocate the arrays required for crash
  crsh_r_k.resize(numRow);
  crsh_c_k.resize(numCol);

  // Add one since need [0..crsh_mx_pri_v]
  crsh_r_pri_k_hdr.resize((crsh_mx_pri_v + 1) * (numCol + 1));
  crsh_r_pri_k_lkf.resize(numRow);
  crsh_r_pri_k_lkb.resize(numRow);
  crsh_r_pri_mn_r_k.resize(crsh_mx_pri_v + 1);

  crsh_mtx_c_mx_abs_v.resize(numCol);
  CrshARvalue.resize(numEl);
  CrshARindex.resize(numEl);
  CrshARstart.resize(numRow + 1);
  crsh_act_r.resize(numRow);
  crsh_act_c.resize(numCol);

  crsh_iz_vr_ty();

  if (crash_strategy == SIMPLEX_CRASH_STRATEGY_BASIC) {
    // For the basis crash, once the row and column priorities have
    // been set, start from a logical basis
    for (int iCol = 0; iCol < numCol; iCol++)
      simplex_basis.nonbasicFlag_[iCol] = NONBASIC_FLAG_TRUE;
    for (int iRow = 0; iRow < numRow; iRow++)
      simplex_basis.nonbasicFlag_[numCol + iRow] = NONBASIC_FLAG_FALSE;
  }
  mx_r_pri = crsh_mn_pri_v;
  for (int r_n = 0; r_n < numRow; r_n++) {
    mx_r_pri = max(mx_r_pri, crsh_r_ty_pri_v[crsh_r_ty[r_n]]);
  }
  mx_c_pri = crsh_mn_pri_v;
  for (int c_n = 0; c_n < numCol; c_n++) {
    mx_c_pri = max(mx_c_pri, crsh_c_ty_pri_v[crsh_c_ty[c_n]]);
  }

  if ((!alw_al_bs_cg) && (mx_r_pri + mx_c_pri <= crsh_mx_pri_v)) return;
  for (int c_n = 0; c_n < numCol + 1; c_n++) {
    for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
      crsh_r_pri_k_hdr[pri_v * (numCol + 1) + c_n] = no_lk;
    }
  }
  mn_abs_pv_v = HIGHS_CONST_INF;
  mn_rlv_pv_v = HIGHS_CONST_INF;
  n_abs_pv_no_ok = 0;
  n_rlv_pv_no_ok = 0;
  // Determine the status and type of each row
  for (int r_n = 0; r_n < numRow; r_n++) {
    int vr_ty = crsh_r_ty[r_n];
    int pri_v = crsh_r_ty_pri_v[vr_ty];
    if (pri_v == crsh_no_act_pri_v) {
      // Rows with no priority value are free - and should be removed
      // by presolve. They will be basic
      crsh_act_r[r_n] = crsh_vr_st_no_act;
    } else {
      crsh_act_r[r_n] = crsh_vr_st_act;
    }
    // Initialise the count of active entries in the row
    crsh_r_k[r_n] = 0;
  }

  // Determine the status and type of each column and compute the
  // number of active entries in each row and column
  for (int c_n = 0; c_n < numCol; c_n++) {
    int vr_ty = crsh_c_ty[c_n];
    int pri_v = crsh_c_ty_pri_v[vr_ty];
    int c_n_en = Astart[c_n + 1] - Astart[c_n];
    if ((pri_v == crsh_no_act_pri_v) || (c_n_en == 0)) {
      // Columns with no priority value are fixed or zero - and should
      // be removed by presolve. They will be nonbasic
      crsh_act_c[c_n] = crsh_vr_st_no_act;
    } else {
      crsh_act_c[c_n] = crsh_vr_st_act;
    }
    //  TODO Count the original number of columns of this type
    if (crsh_act_c[c_n] == crsh_vr_st_act) {
      // The column is active: count the active rows and update the
      // number of active entries in those rows
      //
      // Initialise the count of active entries in the row
      crsh_c_k[c_n] = 0;
      for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
        int r_n = Aindex[el_n];
        if (crsh_act_r[r_n] != crsh_vr_st_no_act) {
          // The row is active so increase the number of active
          // entries in the column and row
          crsh_c_k[c_n] += 1;
          crsh_r_k[r_n] += 1;
        }
      }
    }
  }
  // Now that the row counts are known, make any zero rows non-active.
  //
  // Form linked list of active rows with each priority.
  //
  for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
    crsh_r_pri_mn_r_k[pri_v] = numCol + 1;
  }

  for (int r_n = 0; r_n < numRow; r_n++) {
    if (crsh_act_r[r_n] == crsh_vr_st_no_act) continue;
    // Row appears active but check whether it is empty
    if (crsh_r_k[r_n] <= 0) {
      crsh_act_r[r_n] = crsh_vr_st_no_act;
      continue;
    }
    // Add as the header of the appropriate priority value list
    int pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
    int r_k = crsh_r_k[r_n];
    // Add as the header of rows of a given priority with this count
    int nx_r_n = crsh_r_pri_k_hdr[pri_v * (numCol + 1) + r_k];
    crsh_r_pri_k_hdr[pri_v * (numCol + 1) + r_k] = r_n;
    crsh_r_pri_k_lkb[r_n] = no_lk;
    crsh_r_pri_k_lkf[r_n] = nx_r_n;
    if (nx_r_n != no_lk) crsh_r_pri_k_lkb[nx_r_n] = r_n;
    // Update the minimum row count for this priority value
    crsh_r_pri_mn_r_k[pri_v] = min(r_k, crsh_r_pri_mn_r_k[pri_v]);
  }
  // Set up the row-wise matrix. Although not essential, it is
  // convenient to avoid columns which are not active since the row
  // counts determined above can be used.
  CrshARstart[0] = 0;
  for (int r_n = 0; r_n < numRow; r_n++) {
    CrshARstart[r_n + 1] = CrshARstart[r_n] + crsh_r_k[r_n];
  }
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_act_c[c_n] == crsh_vr_st_no_act) continue;
    // Column is active so find its largest entry in updating the row counts
    crsh_mtx_c_mx_abs_v[c_n] = 0.0;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      crsh_mtx_c_mx_abs_v[c_n] =
          max(fabs(Avalue[el_n]), crsh_mtx_c_mx_abs_v[c_n]);
      if (crsh_act_r[r_n] == crsh_vr_st_no_act) continue;
      int r_el_n = CrshARstart[r_n];
      CrshARindex[r_el_n] = c_n;
      CrshARvalue[r_el_n] = Avalue[el_n];
      CrshARstart[r_n] += 1;
    }
  }
  for (int r_n = 0; r_n < numRow; r_n++) {
    CrshARstart[r_n] -= crsh_r_k[r_n];
  }
}

#ifdef HiGHSDEV
void HCrash::ltssf_ck_da() {
  bool er_fd = false;
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (crsh_act_r[r_n] == crsh_vr_st_no_act) continue;
    int k = crsh_r_k[r_n];
    int ck_k = 0;
    for (int el_n = CrshARstart[r_n]; el_n < CrshARstart[r_n + 1]; el_n++) {
      int c_n = CrshARindex[el_n];
      if (crsh_act_c[c_n] != crsh_vr_st_no_act) ck_k += 1;
    }
    if (ck_k != k) {
      er_fd = true;
      printf(
          "ERROR: Row %d has number of entries error: True = %d; Updated = "
          "%d\n",
          r_n, ck_k, k);
    }
  }
  for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
    int mn_r_k = numCol + 1;
    for (int k = 1; k < numCol + 1; k++) {
      // Check the rows with this priority value and count
      int r_n = crsh_r_pri_k_hdr[pri_v * (numCol + 1) + k];
      if (r_n == no_ix) continue;
      do {
        int ck_pri_v = crsh_r_ty_pri_v[crsh_r_ty[r_n]];
        if (ck_pri_v != pri_v) {
          er_fd = true;
          printf("ERROR: Row %d has ck_pri_v = %d but pri_v = %d\n", r_n,
                 ck_pri_v, pri_v);
        }
        int ck_k = crsh_r_k[r_n];
        if (ck_k != k) {
          er_fd = true;
          printf("ERROR: Row %d has ck_k = %d but k = %d\n", r_n, ck_k, k);
        }
        int nx_r_n = crsh_r_pri_k_lkf[r_n];
        if (nx_r_n != no_lk) {
          int prev_nx_r_n = crsh_r_pri_k_lkb[nx_r_n];
          if (prev_nx_r_n != r_n) {
            er_fd = true;
            printf(
                "ERROR: Back link error for nx_r_n = %d: prev_nx_r_n = %d but "
                "r_n = %d\n",
                nx_r_n, prev_nx_r_n, r_n);
          }
        }
        // Update the true minimum row count
        mn_r_k = min(k, mn_r_k);
        r_n = nx_r_n;
      } while (r_n != no_lk);
    }
    if (crsh_r_pri_mn_r_k[pri_v] != mn_r_k)
      printf(
          "ERROR: Priority %d has crsh_r_pri_mn_r_k = %d < %d = true mn_r_k\n",
          pri_v, crsh_r_pri_mn_r_k[pri_v], mn_r_k);
    crsh_r_pri_mn_r_k[pri_v] = mn_r_k;
  }
  if (er_fd) {
    printf("Error in LTSSF data\n");
  } else {
    printf("No error in LTSSF data\n");
  }
}
#endif

void HCrash::ltssf_cz_r() {
  cz_r_n = no_ix;
  //  Choose a row from the active submatrix
  if (crsh_fn_cf_pri_v > crsh_fn_cf_k) {
    // Search according to maximum priority function
    for (int pri_v = crsh_mx_pri_v; pri_v > crsh_mn_pri_v; pri_v--) {
      int r_k = crsh_r_pri_mn_r_k[pri_v];
      if (r_k > numCol) continue;
      cz_r_n = crsh_r_pri_k_hdr[pri_v * (numCol + 1) + r_k];
      if (cz_r_n == no_ix) {
        printf(
            "ERROR: header for pri_v = %d and count = %d is empty for "
            "crsh_r_pri_mn_r_k[pri_v] = %d\n",
            pri_v, r_k, crsh_r_pri_mn_r_k[pri_v]);
      }
      break;
    }
  } else {
    // Search according to number of active entries
    int mn_r_k = numCol + 1;
    for (int pri_v = crsh_mx_pri_v; pri_v > crsh_mn_pri_v; pri_v--) {
      int r_k = crsh_r_pri_mn_r_k[pri_v];
      if (r_k < mn_r_k) {
        cz_r_n = crsh_r_pri_k_hdr[pri_v * (numCol + 1) + r_k];
        if (cz_r_n == no_ix) {
          printf(
              "ERROR: header for pri_v = %d and count = %d is empty for "
              "crsh_r_pri_mn_r_k[pri_v] = %d\n",
              pri_v, r_k, crsh_r_pri_mn_r_k[pri_v]);
        }
        mn_r_k = r_k;
        if (mn_r_k == 1) break;
      }
    }
  }
}

void HCrash::ltssf_cz_c() {
  HighsLp& simplex_lp = workHMO.simplex_lp_;
  const int objSense = (int)simplex_lp.sense_;
  const double* colCost = &simplex_lp.colCost_[0];

  cz_c_n = no_ix;
  int su_r_c_pri_v_lm = crsh_mx_pri_v;
  // Choose a column which has maximium priority function amongst those
  // with entries in the selected row---making sure that the pivot is
  // acceptable numerically.
  //
  // Note
  //
  // Numerical checking is switched off by setting the tolerances to
  // zero.
  //
  // Make sure that tl_crsh_rlv_pv_v < 1 otherwise test
  // abs_c_v/crsh_mtx_c_mx_abs_v(c_n) .gt. tl_crsh_rlv_pv_v will never
  // be satisfied

  if (alw_al_bs_cg)
    su_r_c_pri_v_lm = -crsh_mx_pri_v;
  else
    su_r_c_pri_v_lm = crsh_mx_pri_v;

  n_eqv_c = 0;
  pv_v = 0.0;
  double mn_co = HIGHS_CONST_INF;
  int mx_c_pri_fn_v = -HIGHS_CONST_I_INF;
  for (int el_n = CrshARstart[cz_r_n]; el_n < CrshARstart[cz_r_n + 1]; el_n++) {
    int c_n = CrshARindex[el_n];
    if (crsh_act_c[c_n] == crsh_vr_st_no_act) continue;
    // Don't allow the row to be replaced by a column whose priority
    // to remain nonbasic is the same or greater.
    if (!alw_al_bs_cg &&
        (crsh_c_ty_pri_v[crsh_c_ty[c_n]] + cz_r_pri_v <= su_r_c_pri_v_lm))
      continue;
    // If column is worse than current best then break
    int c_pri_fn_v = crsh_fn_cf_pri_v * crsh_c_ty_pri_v[crsh_c_ty[c_n]] -
                     crsh_fn_cf_k * crsh_c_k[c_n];
    if (c_pri_fn_v < mx_c_pri_fn_v) continue;
    // Pivot is OK if pivots are not being checked
    bool pv_ok = no_ck_pv;
    if (!no_ck_pv) {
      // Check the matrix entry if it may be used as a pivot.
      double abs_c_v = fabs(CrshARvalue[el_n]);
      bool abs_pv_v_ok = abs_c_v > tl_crsh_abs_pv_v;
      bool rlv_pv_v_ok = abs_c_v > tl_crsh_rlv_pv_v * crsh_mtx_c_mx_abs_v[c_n];
      if (!abs_pv_v_ok) n_abs_pv_no_ok += 1;
      if (!rlv_pv_v_ok) n_rlv_pv_no_ok += 1;
      pv_ok = abs_pv_v_ok && rlv_pv_v_ok;
    }
    // If the pivot is not OK then break
    if (!pv_ok) continue;
    // Pivot is OK [or not being checked]
    if (c_pri_fn_v > mx_c_pri_fn_v) {
      // Better column than best so far
      mx_c_pri_fn_v = c_pri_fn_v;
      cz_c_n = c_n;
      pv_v = CrshARvalue[el_n];
      n_eqv_c = 1;
      double sense_col_cost = objSense * colCost[c_n];
      mn_co = sense_col_cost;
    } else if (c_pri_fn_v == mx_c_pri_fn_v) {
      // Matches the best column so far: possibly break tie on cost
      double sense_col_cost = objSense * colCost[c_n];
      if (mn_co_tie_bk && (sense_col_cost < mn_co)) {
        cz_c_n = c_n;
        pv_v = CrshARvalue[el_n];
        n_eqv_c = 1;
        mn_co = sense_col_cost;
      }
      n_eqv_c += 1;
    }
  }
}

#ifdef HiGHSDEV
void HCrash::tsSing() {
  printf("\nTesting singularity Crash\n");
  int nBcVr = 0;
  // Make columns basic until they are either all basic or the number
  // of basic variables reaches numRow
  for (int c_n = 0; c_n < numTot; c_n++) {
    int r_n = c_n;
    int columnIn = c_n;
    int columnOut = numCol + r_n;
    workHMO.simplex_basis_.nonbasicFlag_[columnIn] = NONBASIC_FLAG_FALSE;
    workHMO.simplex_basis_.nonbasicFlag_[columnOut] = NONBASIC_FLAG_TRUE;
    nBcVr++;
    if (nBcVr == numRow) break;
  }
}

void HCrash::ltssf_rp_r_k() {
  printf("\n        : ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", r_n);
  }
  printf(" \n");
  printf("k: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_k[r_n]);
  }
  printf(" \n");
}

void HCrash::ltssf_rp_r_pri() {
  printf("\nRow       : ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", r_n);
  }
  printf(" \n");
  printf("crsh_r_ty:  ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_ty[r_n]);
  }
  printf(" \n");
  printf("pri_v:      ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_ty_pri_v[crsh_r_ty[r_n]]);
  }
  printf(" \n");
  printf("act_r:      ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_act_r[r_n]);
  }
}

void HCrash::ltssf_rp_pri_k_da() {
  for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
    printf("Pri = %1d: ", pri_v);
    for (int k = 1; k < numCol + 1; k++) {
      int hdr_ix = pri_v * (numCol + 1) + k;
      printf(" %2d", crsh_r_pri_k_hdr[hdr_ix]);
    }
    printf("\n");
  }
  printf("r_n: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", r_n);
  }
  printf("\n");
  printf("Act: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_act_r[r_n]);
  }
  printf("\n");
  printf("Cnt: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_k[r_n]);
  }
  printf("\n");
  printf("Lkf: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_pri_k_lkf[r_n]);
  }
  printf("\n");
  printf("Lkb: ");
  for (int r_n = 0; r_n < numRow; r_n++) {
    printf(" %2d", crsh_r_pri_k_lkb[r_n]);
  }
  printf("\n");
  for (int pri_v = crsh_mn_pri_v; pri_v < crsh_mx_pri_v + 1; pri_v++) {
    for (int k = 1; k < numCol + 1; k++) {
      int hdr_ix = pri_v * (numCol + 1) + k;
      int r_n = crsh_r_pri_k_hdr[hdr_ix];
      //   printf("Pri = %1d; k = %1d: hdr = %d\n", pri_v,
      // k, r_n);
      if (r_n == no_lk) continue;
      int nx_r_n = crsh_r_pri_k_lkf[r_n];
      printf("Pri = %1d; k = %1d: %3d %3d <-  ", pri_v, k, r_n, nx_r_n);
      r_n = nx_r_n;
      if (r_n != no_lk) {
        do {
          int prev_r_n = crsh_r_pri_k_lkb[r_n];
          int nx_r_n = crsh_r_pri_k_lkf[r_n];
          printf("  <-%3d %3d %3d ->  ", prev_r_n, r_n, nx_r_n);
          r_n = nx_r_n;
        } while (r_n != no_lk);
      }
      printf("\n");
    }
  }
}

#endif

void HCrash::crsh_iz_vr_ty() {
  HighsLp& simplex_lp = workHMO.simplex_lp_;
  const double* colLower = &simplex_lp.colLower_[0];
  const double* colUpper = &simplex_lp.colUpper_[0];
  const double* rowLower = &simplex_lp.rowLower_[0];
  const double* rowUpper = &simplex_lp.rowUpper_[0];
  const int* nonbasicFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  // Allocate the arrays required for crash
  crsh_r_ty.resize(numRow);
  crsh_c_ty.resize(numCol);
  if (crash_strategy == SIMPLEX_CRASH_STRATEGY_BASIC) {
    for (int r_n = 0; r_n < numRow; r_n++) {
      if (nonbasicFlag[numCol + r_n] == NONBASIC_FLAG_TRUE)
        crsh_r_ty[r_n] = crsh_vr_ty_non_bc;
      else
        crsh_r_ty[r_n] = crsh_vr_ty_bc;
    }
    for (int c_n = 0; c_n < numCol; c_n++) {
      if (nonbasicFlag[c_n] == NONBASIC_FLAG_TRUE)
        crsh_c_ty[c_n] = crsh_vr_ty_non_bc;
      else
        crsh_c_ty[c_n] = crsh_vr_ty_bc;
    }
  } else {
    for (int r_n = 0; r_n < numRow; r_n++) {
      if (rowUpper[r_n] >= HIGHS_CONST_INF) {
        if (rowLower[r_n] <= -HIGHS_CONST_INF)
          crsh_r_ty[r_n] = crsh_vr_ty_fr;  // Free row
        else
          crsh_r_ty[r_n] = crsh_vr_ty_1_sd;  // Lower-bounded (1-sided) row
      } else {
        if (rowLower[r_n] <= -HIGHS_CONST_INF)
          crsh_r_ty[r_n] = crsh_vr_ty_1_sd;  // Upper-bonded (1-sided) row
        else {
          // Two-sided row - maybe fixed (equality)
          if (rowLower[r_n] != rowUpper[r_n])
            crsh_r_ty[r_n] = crsh_vr_ty_2_sd;  // 2-sided row
          else
            crsh_r_ty[r_n] = crsh_vr_ty_fx;  // Fixed (equality) row
        }
      }
    }
    // Set up the column variable types for crash
    for (int c_n = 0; c_n < numCol; c_n++) {
      if (colUpper[c_n] >= HIGHS_CONST_INF) {
        if (colLower[c_n] <= -HIGHS_CONST_INF)
          crsh_c_ty[c_n] = crsh_vr_ty_fr;  // Free column
        else
          crsh_c_ty[c_n] = crsh_vr_ty_1_sd;  // Lower-bounded (1-sided) column
      } else {
        if (colLower[c_n] <= -HIGHS_CONST_INF)
          crsh_c_ty[c_n] = crsh_vr_ty_1_sd;  // Upper-bonded (1-sided) column
        else {
          // Two-sided row - maybe fixed (equality)
          if (colLower[c_n] != colUpper[c_n])
            crsh_c_ty[c_n] = crsh_vr_ty_2_sd;  // 2-sided column
          else
            crsh_c_ty[c_n] = crsh_vr_ty_fx;  // Fixed column
        }
      }
    }
  }
#ifdef HiGHSDEV
  // Allocate the arrays to analyse crash
  crsh_vr_ty_og_n_r.resize(crsh_num_vr_ty);
  crsh_vr_ty_rm_n_r.resize(crsh_num_vr_ty);
  crsh_vr_ty_og_n_c.resize(crsh_num_vr_ty);
  crsh_vr_ty_add_n_c.resize(crsh_num_vr_ty);
  crsh_bs_vr_ty_n_r.resize(crsh_num_vr_ty);
  crsh_bs_vr_ty_n_c.resize(crsh_num_vr_ty);
  crsh_nonbc_vr_ty_n_r.resize(crsh_num_vr_ty);
  crsh_nonbc_vr_ty_n_c.resize(crsh_num_vr_ty);
  // Initialise the counts of numbers and changes of variable types -
  // just for reporting
  for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_num_vr_ty; vr_ty++) {
    crsh_vr_ty_og_n_r[vr_ty] = 0;
    crsh_vr_ty_og_n_c[vr_ty] = 0;
    crsh_vr_ty_rm_n_r[vr_ty] = 0;
    crsh_vr_ty_add_n_c[vr_ty] = 0;
    crsh_bs_vr_ty_n_r[vr_ty] = 0;
    crsh_bs_vr_ty_n_c[vr_ty] = 0;
    crsh_nonbc_vr_ty_n_r[vr_ty] = 0;
    crsh_nonbc_vr_ty_n_c[vr_ty] = 0;
  }
  for (int r_n = 0; r_n < numRow; r_n++) crsh_vr_ty_og_n_r[crsh_r_ty[r_n]] += 1;
  for (int c_n = 0; c_n < numCol; c_n++) crsh_vr_ty_og_n_c[crsh_c_ty[c_n]] += 1;
#endif
}

#ifdef HiGHSDEV
void HCrash::crsh_an_c_co() {
  HighsLp& simplex_lp = workHMO.simplex_lp_;
  const int objSense = (int)simplex_lp.sense_;
  const double* colCost = &simplex_lp.colCost_[0];
  const double* colLower = &simplex_lp.colLower_[0];
  const double* colUpper = &simplex_lp.colUpper_[0];

  int n_ze_c_co = 0;
  int n_fs_c_co = 0;

  for (int c_n = 0; c_n < numCol; c_n++) {
    double sense_col_cost = objSense * colCost[c_n];
    if (sense_col_cost == 0.0) {
      n_ze_c_co += 1;
      n_fs_c_co += 1;
      continue;
    }
    if (colUpper[c_n] >= HIGHS_CONST_INF) {
      // Free column: nonzero cost cannot be feasible
      if (colLower[c_n] > -HIGHS_CONST_INF) {
        // Lower-bounded (1-sided) column: non-negative cost is feasible
        double sense_col_cost = objSense * colCost[c_n];
        if (sense_col_cost >= 0.0) n_fs_c_co += 1;
      }
    } else {
      if (colLower[c_n] <= -HIGHS_CONST_INF) {
        // Upper-bonded (1-sided) column: non-positive cost is feasible
        double sense_col_cost = objSense * colCost[c_n];
        if (sense_col_cost <= 0.0) n_fs_c_co += 1;
      } else {
        // Two-sided column: any cost is feasible
        n_fs_c_co += 1;
      }
    }
  }
  printf(" Model has %7d Ze costs (%3d%%)\n", n_ze_c_co,
         (100 * n_ze_c_co) / numCol);
  printf(" Model has %7d Fs costs (%3d%%)\n", n_fs_c_co,
         (100 * n_fs_c_co) / numCol);
}

void HCrash::crsh_rp_r_c_st(const int mode) {
  string TyNm;
  int ck_su_n_c = 0;
  int ck_su_n_r = 0;
  int ck_su_n_bc_vr = 0;
  int ck_su_n_nonbc_vr = 0;
  int n_ps = 2;
  if (mode == 1) n_ps = 1;
  for (int ps_n = 0; ps_n < n_ps; ps_n++) {
    if (ps_n == 1) {
      if (mode == 0)
        printf("grep_CharCrash,Rows");
      else if (mode == 2)
        printf("grep_CharCrash,Basic");
      else
        printf("grep_CharCrash,Nonbasic");
      for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_num_vr_ty; vr_ty++) {
        TyNm = crsh_nm_o_crsh_vr_ty(vr_ty);
        if (mode == 0) {
          printf(",%s", TyNm.c_str());
        } else {
          printf(",%s_Row", TyNm.c_str());
          printf(",%s_Col", TyNm.c_str());
        }
      }
      printf("\n");
      if (mode == 0)
        printf("grep_CharCrash,%d", numRow);
      else if (mode == 2)
        printf("grep_CharCrash,%d", numRow);
      else if (mode == 3)
        printf("grep_CharCrash,%d", numCol);
    }
    for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_num_vr_ty; vr_ty++) {
      TyNm = crsh_nm_o_crsh_vr_ty(vr_ty);
      if (mode == 0) {
        if (ps_n == 0) ck_su_n_r += crsh_vr_ty_og_n_r[vr_ty];
        int lc_pct = (100 * crsh_vr_ty_og_n_r[vr_ty]) / numRow;
        if (ps_n == 0) {
          if (crsh_vr_ty_og_n_r[vr_ty] > 0)
            printf(" Model has %7d %3s rows (%3d%%)\n",
                   crsh_vr_ty_og_n_r[vr_ty], TyNm.c_str(), lc_pct);
        } else {
          printf(",%7d", crsh_vr_ty_og_n_r[vr_ty]);
        }
      } else if (mode == 1) {
        if (crsh_vr_ty_og_n_r[vr_ty] > 0)
          printf(" Removed %7d of %7d %3s rows (%3d%%)\n",
                 crsh_vr_ty_rm_n_r[vr_ty], crsh_vr_ty_og_n_r[vr_ty],
                 TyNm.c_str(),
                 (100 * crsh_vr_ty_rm_n_r[vr_ty]) / crsh_vr_ty_og_n_r[vr_ty]);
      } else if (mode == 2) {
        if (ps_n == 0) {
          ck_su_n_bc_vr += crsh_bs_vr_ty_n_r[vr_ty];
          ck_su_n_bc_vr += crsh_bs_vr_ty_n_c[vr_ty];
          ck_su_n_nonbc_vr += crsh_nonbc_vr_ty_n_r[vr_ty];
          ck_su_n_nonbc_vr += crsh_nonbc_vr_ty_n_c[vr_ty];
          if (crsh_bs_vr_ty_n_r[vr_ty] > 0)
            printf(" Basic    variables contain %7d %3s rows (%3d%%)\n",
                   crsh_bs_vr_ty_n_r[vr_ty], TyNm.c_str(),
                   (100 * crsh_bs_vr_ty_n_r[vr_ty]) / numRow);
          if (crsh_bs_vr_ty_n_c[vr_ty] > 0)
            printf(" Basic    variables contain %7d %3s cols (%3d%%)\n",
                   crsh_bs_vr_ty_n_c[vr_ty], TyNm.c_str(),
                   (100 * crsh_bs_vr_ty_n_c[vr_ty]) / numRow);
        } else {
          printf(",%d,%d", crsh_bs_vr_ty_n_r[vr_ty], crsh_bs_vr_ty_n_c[vr_ty]);
        }
      } else {
        if (ps_n == 0) {
          if (crsh_nonbc_vr_ty_n_c[vr_ty] > 0)
            printf(" Nonbasic variables contain %7d %3s cols (%3d%%)\n",
                   crsh_nonbc_vr_ty_n_c[vr_ty], TyNm.c_str(),
                   (100 * crsh_nonbc_vr_ty_n_c[vr_ty]) / numCol);
          if (crsh_nonbc_vr_ty_n_r[vr_ty] > 0)
            printf(" Nonbasic variables contain %7d %3s rows (%3d%%)\n",
                   crsh_nonbc_vr_ty_n_r[vr_ty], TyNm.c_str(),
                   (100 * crsh_nonbc_vr_ty_n_r[vr_ty]) / numCol);
        } else {
          printf(",%d,%d", crsh_nonbc_vr_ty_n_r[vr_ty],
                 crsh_nonbc_vr_ty_n_c[vr_ty]);
        }
      }
    }
    if (ps_n == 1) printf("\n");
  }
  if (mode == 0) assert(ck_su_n_r == numRow);
  if (mode == 2) assert(ck_su_n_bc_vr == numRow);
  if (mode == 2) assert(ck_su_n_nonbc_vr == numCol);
  if (mode <= 1) {
    for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_num_vr_ty; vr_ty++) {
      TyNm = crsh_nm_o_crsh_vr_ty(vr_ty);
      if (mode == 0) ck_su_n_c += crsh_vr_ty_og_n_c[vr_ty];
      if (crsh_vr_ty_og_n_c[vr_ty] > 0)
        printf(" Model has %7d %3s cols (%3d%%)\n", crsh_vr_ty_og_n_c[vr_ty],
               TyNm.c_str(), (100 * crsh_vr_ty_og_n_c[vr_ty]) / numCol);

      else if (mode == 1) {
        if (crsh_vr_ty_og_n_c[vr_ty] > 0)
          printf(" Added   %7d of %7d %3s cols (%3d%%)\n",
                 crsh_vr_ty_add_n_c[vr_ty], crsh_vr_ty_og_n_c[vr_ty],
                 TyNm.c_str(),
                 (100 * crsh_vr_ty_add_n_c[vr_ty]) / crsh_vr_ty_og_n_c[vr_ty]);
      }
    }
    if (mode == 0) assert(ck_su_n_c == numCol);
  }
}
void HCrash::crsh_an_r_c_st_af() {
  const int* Astart = &workHMO.simplex_lp_.Astart_[0];
  for (int k = 0; k < numRow; k++) {
    int vr_n = workHMO.simplex_basis_.basicIndex_[k];
    if (vr_n < numCol) {
      int c_n = vr_n;
      crsh_bs_vr_ty_n_c[crsh_c_ty[c_n]] += 1;
    } else {
      int r_n = vr_n - numCol;
      crsh_bs_vr_ty_n_r[crsh_r_ty[r_n]] += 1;
    }
  }

  for (int vr_n = 0; vr_n < numTot; vr_n++) {
    if (workHMO.simplex_basis_.nonbasicFlag_[vr_n] == 0) continue;
    if (vr_n < numCol) {
      int c_n = vr_n;
      crsh_nonbc_vr_ty_n_c[crsh_c_ty[c_n]] += 1;
    } else {
      int r_n = vr_n - numCol;
      crsh_nonbc_vr_ty_n_r[crsh_r_ty[r_n]] += 1;
    }
  }
  int bs_mtx_n_struc_el = 0;
  for (int r_n = 0; r_n < numRow; r_n++) {
    int vr_n = workHMO.simplex_basis_.basicIndex_[r_n];
    if (vr_n < numCol) {
      int c_n_el = Astart[vr_n + 1] - Astart[vr_n];
      bs_mtx_n_struc_el += c_n_el;
    }
  }

  crsh_rp_r_c_st(1);
  crsh_rp_r_c_st(2);
  printf(" Basis    matrix    contains%7d structural entries\n",
         bs_mtx_n_struc_el);
  crsh_rp_r_c_st(3);
}

string HCrash::crsh_nm_o_crsh_vr_ty(const int vr_ty) {
  string TyNm;
  if (crash_strategy == SIMPLEX_CRASH_STRATEGY_BASIC) {
    if (vr_ty == crsh_vr_ty_non_bc)
      TyNm = "NBc";
    else if (vr_ty == crsh_vr_ty_bc)
      TyNm = " Bc";
    else
      printf("Unrecognised type %d\n", vr_ty);
  } else {
    if (vr_ty == crsh_vr_ty_fx)
      TyNm = "Fx ";
    else if (vr_ty == crsh_vr_ty_2_sd)
      TyNm = "2sd";
    else if (vr_ty == crsh_vr_ty_1_sd)
      TyNm = "1sd";
    else if (vr_ty == crsh_vr_ty_fr)
      TyNm = "Fr ";
    else
      printf("Unrecognised type %d\n", vr_ty);
  }
  return TyNm;
}

#endif
