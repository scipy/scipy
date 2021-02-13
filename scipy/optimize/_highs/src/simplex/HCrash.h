/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HCrash.h
 * @brief Bixby and Maros-style crash for the HiGHS simplex solver
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HCRASH_H_
#define SIMPLEX_HCRASH_H_

#include <string>
#include <vector>

#include "lp_data/HighsModelObject.h"

class HMatrix;

// LTSSF scalar parameters
const int crsh_vr_st_no_act = 0;
const int crsh_vr_st_act = 1;

// Crash variable types
// Basis-preserving crash:
const int crsh_vr_ty_non_bc = 0;
const int crsh_vr_ty_bc = 1;
// Standard crash:
const int crsh_vr_ty_fx = 0;
const int crsh_vr_ty_2_sd = 1;
const int crsh_vr_ty_1_sd = 2;
const int crsh_vr_ty_fr = 3;

// Null header for linked lists
const int no_lk = -1;

// Null value for chosen row/column index
const int no_ix = no_lk;

// LTSSF scalar control parameters
const double tl_crsh_abs_pv_v = 1e-4;
const double tl_crsh_rlv_pv_v = 1e-2;
// Switches for LTSSF checking and reporting
const int ltssf_ck_fq = 0;
#ifdef HiGHSDEV
const bool reportCrashData = false;
const bool reportBixbyPass = false;
#endif

/**
 * @brief Bixby and Maros-style crash for the HiGHS simplex solver
 */
class HCrash {
 public:
  HCrash(HighsModelObject& model_object) : workHMO(model_object) {}
  /**
   * @brief Determine a particular crash basis for a given model instance
   */
  void crash(const int pass_crash_strategy);

 private:
  // Internal methods
  void bixby();
  bool bixby_iz_da();
  void bixby_rp_mrt();

  void ltssf();
  void ltssf_iz_mode();
  void ltssf_iz_da();
  void ltssf_iterate();
  void ltssf_u_da();
  void ltssf_u_da_af_bs_cg();
  void ltssf_u_da_af_no_bs_cg();
#ifdef HiGHSDEV
  void ltssf_ck_da();
#endif
  void ltssf_cz_r();
  void ltssf_cz_c();
#ifdef HiGHSDEV
  void tsSing();
  void ltssf_rp_r_k();
  void ltssf_rp_r_pri();
  void ltssf_rp_pri_k_da();
#endif

  void crsh_iz_vr_ty();

#ifdef HiGHSDEV
  void crsh_an_c_co();
  void crsh_rp_r_c_st(const int mode);
  void crsh_an_r_c_st_af();
  std::string crsh_nm_o_crsh_vr_ty(const int vr_ty);
#endif

  // Model to be crashed
  HighsModelObject& workHMO;

  // Crash strategy to be used
  int crash_strategy;

  // Model dimensions
  int numCol;
  int numRow;
  int numTot;

  //    LTSSF arrays
  std::vector<int> crsh_r_ty_pri_v;
  std::vector<int> crsh_c_ty_pri_v;
  std::vector<int> crsh_r_ty;
  std::vector<int> crsh_c_ty;
  std::vector<int> crsh_r_k;
  std::vector<int> crsh_c_k;

  std::vector<int> crsh_r_pri_k_hdr;
  std::vector<int> crsh_r_pri_k_lkf;
  std::vector<int> crsh_r_pri_k_lkb;
  std::vector<int> crsh_r_pri_mn_r_k;

  std::vector<int> crsh_r_pri_hdr;
  std::vector<int> crsh_r_pri_lkb;
  std::vector<int> crsh_r_pri_lkf;

  std::vector<int> crsh_r_k_hdr;
  std::vector<int> crsh_r_k_lkb;
  std::vector<int> crsh_r_k_lkf;

#ifdef HiGHSDEV
  std::vector<int> crsh_vr_ty_og_n_r;
  std::vector<int> crsh_vr_ty_rm_n_r;
  std::vector<int> crsh_vr_ty_og_n_c;
  std::vector<int> crsh_vr_ty_add_n_c;

  std::vector<int> crsh_bs_vr_ty_n_r;
  std::vector<int> crsh_bs_vr_ty_n_c;
  std::vector<int> crsh_nonbc_vr_ty_n_r;
  std::vector<int> crsh_nonbc_vr_ty_n_c;
#endif

  std::vector<double> crsh_mtx_c_mx_abs_v;
  std::vector<double> CrshARvalue;
  std::vector<int> CrshARindex;
  std::vector<int> CrshARstart;
  std::vector<int> crsh_act_r;
  std::vector<int> crsh_act_c;

  std::vector<double> bixby_mrt_v;
  std::vector<double> heap_v;
  std::vector<double> bixby_pseudo_pv_v;
  std::vector<int> bixby_mrt_ix;
  std::vector<int> heap_ix;
  std::vector<int> bixby_pv_in_r;
  std::vector<int> bixby_vr_in_r;
  std::vector<int> bixby_r_k;
  // std::vector<int> bixby_ze_r_k;

  // LTSSF scalar identifiers
  // int crsh_mode;
  int crsh_f_vr_ty;
  int crsh_l_vr_ty;
  int crsh_num_vr_ty;

  int crsh_mn_pri_v;      // = 0;
  int crsh_mx_pri_v;      // = 3;
  int crsh_no_act_pri_v;  // = crsh_mn_pri_v;

  int crsh_fn_cf_pri_v;
  int crsh_fn_cf_k;
  bool mn_co_tie_bk;
  bool alw_al_bs_cg;
  bool no_ck_pv;
  double bixby_mu_a;
  double bixby_mu_b;

  // LTSSF scalar identifiers
  int n_crsh_ps;
  int n_crsh_bs_cg;
  int cz_r_n;
  int cz_r_pri_v;
  int cz_c_n;
  int n_eqv_c;
  double pv_v;
  double mn_abs_pv_v;
  double mn_rlv_pv_v;
  int mx_r_pri_v;
  int n_abs_pv_no_ok;
  int n_rlv_pv_no_ok;
  int mx_r_pri;
  int mx_c_pri;
  int bixby_n_cdd_r;
  bool bixby_no_nz_c_co;
};

#endif /* SIMPLEX_HCRASH_H_ */
