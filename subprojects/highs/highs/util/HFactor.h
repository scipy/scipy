/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HFactor.h
 * @brief Basis matrix factorization, update and solves for HiGHS
 */
#ifndef HFACTOR_H_
#define HFACTOR_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsAnalysis.h"
#include "util/HVector.h"
#include "util/HighsSparseMatrix.h"

// Uses max and min for local in-line functions
using std::max;
// using std::min;
using std::vector;

const HighsInt kBuildKernelReturnTimeout = -1;

struct InvertibleRepresentation {
  // Factor L
  std::vector<HighsInt> l_pivot_index;
  std::vector<HighsInt> l_pivot_lookup;
  std::vector<HighsInt> l_start;
  std::vector<HighsInt> l_index;
  std::vector<double> l_value;
  std::vector<HighsInt> lr_start;
  std::vector<HighsInt> lr_index;
  std::vector<double> lr_value;

  // Factor U
  std::vector<HighsInt> u_pivot_lookup;
  std::vector<HighsInt> u_pivot_index;
  std::vector<double> u_pivot_value;

  //  HighsInt u_total_x;
  std::vector<HighsInt> u_start;
  std::vector<HighsInt> u_last_p;
  std::vector<HighsInt> u_index;
  std::vector<double> u_value;

  std::vector<HighsInt> ur_start;
  std::vector<HighsInt> ur_lastp;
  std::vector<HighsInt> ur_space;
  std::vector<HighsInt> ur_index;
  std::vector<double> ur_value;
  std::vector<HighsInt> pf_start;
  std::vector<HighsInt> pf_index;
  std::vector<double> pf_value;
  std::vector<HighsInt> pf_pivot_index;
  std::vector<double> pf_pivot_value;
  void clear();
};

/**
 * @brief Basis matrix factorization, update and solves for HiGHS
 *
 * Class for the following
 *
 * Basis matrix factorization \f$PBQ=LU\f$
 *
 * Update according to \f$B'=B+(\mathbf{a}_q-B\mathbf{e}_p)\mathbf{e}_p^T\f$
 *
 * Solves \f$B\mathbf{x}=\mathbf{b}\f$ (FTRAN) and
 * \f$B^T\mathbf{x}=\mathbf{b}\f$ (BTRAN)
 *
 * HFactor is initialised using HFactor::setup, which takes copies of
 * the pointers to the constraint matrix starts, indices, values and
 * basic column indices.
 *
 * Forming \f$PBQ=LU\f$ (INVERT) is performed using HFactor::build
 *
 * Solving \f$B\mathbf{x}=\mathbf{b}\f$ (FTRAN) is performed using
 * HFactor::ftran
 *
 * Solving \f$B^T\mathbf{x}=\mathbf{b}\f$ (BTRAN) is performed using
 * HFactor::btran
 *
 * Updating the invertible representation of the basis matrix
 * according to \f$B'=B+(\mathbf{a}_q-B\mathbf{e}_p)\mathbf{e}_p^T\f$
 * is performed by HFactor::update. UPDATE requires vectors
 * \f$B^{-1}\mathbf{a}_q\f$ and \f$B^{-T}\mathbf{e}_q\f$, together
 * with the index of the pivotal row.
 *
 * HFactor assumes that the basic column indices are kept up-to-date
 * externally as basis changes take place. INVERT permutes the basic
 * column indices, since these define the order of the solution values
 * after FTRAN, and the assumed order of the RHS before BTRAN
 *
 */
class HFactor {
 public:
  HFactor()
      : build_realTick(0.0),
        build_synthetic_tick(0.0),
        rank_deficiency(0),
        basis_matrix_num_el(0),
        invert_num_el(0),
        kernel_dim(0),
        kernel_num_el(0),
        num_row(0),
        num_col(0),
        num_basic(0),
        a_matrix_valid(false),
        a_start(nullptr),
        a_index(nullptr),
        a_value(nullptr),
        basic_index(nullptr),
        pivot_threshold(0.0),
        pivot_tolerance(0.0),
        highs_debug_level(0),
        time_limit_(0.0),
        use_original_HFactor_logic(false),
        debug_report_(false),
        basis_matrix_limit_size(0),
        update_method(0),
        build_timer_(nullptr),
        nwork(0),
        u_merit_x(0),
        // clang-format off
        u_total_x(0) {};
  // clang-format on

  /**
   * @brief Copy problem size and pointers of constraint matrix, and set
   * up space for INVERT
   *
   * Copy problem size and pointers to coefficient matrix, allocate
   * working buffer for INVERT, allocate space for basis matrix, L, U
   * factor and Update buffer, allocated space for Markowitz matrices,
   * count-link-list, L factor and U factor
   */

  void setup(const HighsSparseMatrix& a_matrix,
             std::vector<HighsInt>& basic_index,
             const double pivot_threshold = kDefaultPivotThreshold,
             const double pivot_tolerance = kDefaultPivotTolerance,
             const HighsInt highs_debug_level = kHighsDebugLevelMin,
             const HighsLogOptions* log_options = NULL);

  void setupGeneral(const HighsSparseMatrix* a_matrix, HighsInt num_basic,
                    HighsInt* basic_index,
                    const double pivot_threshold = kDefaultPivotThreshold,
                    const double pivot_tolerance = kDefaultPivotTolerance,
                    const HighsInt highs_debug_level = kHighsDebugLevelMin,
                    const HighsLogOptions* log_options = NULL);

  void setup(const HighsInt num_col,   //!< Number of columns
             const HighsInt num_row,   //!< Number of rows
             const HighsInt* a_start,  //!< Column starts of constraint matrix
             const HighsInt* a_index,  //!< Row indices of constraint matrix
             const double* a_value,    //!< Row values of constraint matrix
             HighsInt* basic_index,    //!< Indices of basic variables
             const double pivot_threshold =
                 kDefaultPivotThreshold,  //!< Pivoting threshold
             const double pivot_tolerance =
                 kDefaultPivotTolerance,  //!< Min absolute pivot
             const HighsInt highs_debug_level = kHighsDebugLevelMin,
             const HighsLogOptions* log_options = NULL,
             const bool use_original_HFactor_logic = true,
             const HighsInt update_method = kUpdateMethodFt);

  void setupGeneral(
      const HighsInt num_col,    //!< Number of columns
      const HighsInt num_row,    //!< Number of rows
      const HighsInt num_basic,  //!< Number of indices in basic_index
      const HighsInt* a_start,   //!< Column starts of constraint matrix
      const HighsInt* a_index,   //!< Row indices of constraint matrix
      const double* a_value,     //!< Row values of constraint matrix
      HighsInt* basic_index,     //!< Indices of "basic" variables
      const double pivot_threshold =
          kDefaultPivotThreshold,  //!< Pivoting threshold
      const double pivot_tolerance =
          kDefaultPivotTolerance,  //!< Min absolute pivot
      const HighsInt highs_debug_level = kHighsDebugLevelMin,
      const HighsLogOptions* log_options = NULL,
      const bool use_original_HFactor_logic = true,
      const HighsInt update_method = kUpdateMethodFt);

  void setupMatrix(
      const HighsInt* a_start,  //!< Column starts of constraint matrix
      const HighsInt* a_index,  //!< Row indices of constraint matrix
      const double* a_value);   //!< Row values of constraint matrix
  void setupMatrix(const HighsSparseMatrix* a_matrix);
  /**
   * @brief Form \f$PBQ=LU\f$ for basis matrix \f$B\f$ or report degree of rank
   * deficiency.
   *
   * @return 0 if successful, otherwise rank_deficiency>0
   *
   */
  HighsInt build(HighsTimerClock* factor_timer_clock_pointer = NULL);

  /**
   * @brief Solve \f$B\mathbf{x}=\mathbf{b}\f$ (FTRAN)
   */
  void ftranCall(
      HVector& vector,                //!< RHS vector \f$\mathbf{b}\f$
      const double expected_density,  //!< Expected density of the result
      HighsTimerClock* factor_timer_clock_pointer = NULL) const;

  void ftranCall(std::vector<double>& vector,
                 HighsTimerClock* factor_timer_clock_pointer = NULL);

  /**
   * @brief Solve \f$B^T\mathbf{x}=\mathbf{b}\f$ (BTRAN)
   */
  void btranCall(
      HVector& vector,                //!< RHS vector \f$\mathbf{b}\f$
      const double expected_density,  //!< Expected density of the result
      HighsTimerClock* factor_timer_clock_pointer = NULL) const;

  void btranCall(std::vector<double>& vector,
                 HighsTimerClock* factor_timer_clock_pointer = NULL);

  /**
   * @brief Update according to
   * \f$B'=B+(\mathbf{a}_q-B\mathbf{e}_p)\mathbf{e}_p^T\f$
   */
  void update(HVector* aq,     //!< Vector \f$B^{-1}\mathbf{a}_q\f$
              HVector* ep,     //!< Vector \f$B^{-T}\mathbf{e}_p\f$
              HighsInt* iRow,  //!< Index of pivotal row
              HighsInt* hint   //!< Reinversion status
  );

  /**
   * @brief Sets pivoting threshold
   */
  bool setPivotThreshold(
      const double new_pivot_threshold = kDefaultPivotThreshold);

  /**
   * @brief Sets a time limit
   */
  void setTimeLimit(const double time_limit);

  /**
   * @brief Updates instance with respect to new columns in the
   * constraint matrix (assuming columns are nonbasic)
   */
  void addCols(const HighsInt num_new_col);

  /**
   * @brief Updates instance with respect to nonbasic columns in the
   * constraint matrix being deleted
   */
  void deleteNonbasicCols(const HighsInt num_deleted_col);

  /**
   * @brief Updates instance with respect to new rows in the
   * constraint matrix (assuming slacks are basic)
   */
  void addRows(const HighsSparseMatrix* ar_matrix);

  /**
   * @brief Wall clock time for INVERT
   */
  double build_realTick;

  /**
   * @brief The synthetic clock for INVERT
   */
  double build_synthetic_tick;

  // Rank deficiency information

  /**
   * @brief Degree of rank deficiency in \f$B\f$
   */
  HighsInt rank_deficiency;

  /**
   * @brief Rows not pivoted on
   */
  vector<HighsInt> row_with_no_pivot;

  /**
   * @brief (Basis matrix) columns not pivoted on
   */
  vector<HighsInt> col_with_no_pivot;

  /**
   * @brief Variables not pivoted on
   */
  vector<HighsInt> var_with_no_pivot;

  /**
   * @brief Gets basic_index since it is private
   */
  const HighsInt* getBaseIndex() const { return basic_index; }

  /**
   * @brief Gets a_start since it is private
   */
  const HighsInt* getAstart() const { return a_start; }

  /**
   * @brief Gets a_index since it is private
   */
  const HighsInt* getAindex() const { return a_index; }

  /**
   * @brief Gets a_value since it is private
   */
  const double* getAvalue() const { return a_value; }

  void reportLu(const HighsInt l_u_or_both = kReportLuBoth,
                const bool full = true) const;
  void reportAsm() const;

  InvertibleRepresentation getInvert() const;
  void setInvert(const InvertibleRepresentation& invert);

  void setDebugReport(const bool debug_report) {
    this->debug_report_ = debug_report;
  }

  // Information required to perform refactorization of the current
  // basis
  RefactorInfo refactor_info_;

  // Properties of data held in HFactor.h
  HighsInt basis_matrix_num_el;
  HighsInt invert_num_el;
  HighsInt kernel_dim;
  HighsInt kernel_num_el;
  HighsInt num_active_nz_;
  double min_time_bound_;
  double max_time_bound_;

  /**
   * Data of the factor
   */

  // private:
  // Problem size, coefficient matrix and update method
  HighsInt num_row;
  HighsInt num_col;
  HighsInt num_basic;
  double inv_num_row;  // 1.0/num_row

 private:
  bool a_matrix_valid;
  const HighsInt* a_start;
  const HighsInt* a_index;
  const double* a_value;
  HighsInt* basic_index;
  double pivot_threshold;
  double pivot_tolerance;
  HighsInt highs_debug_level;
  double time_limit_;

  struct LogData {
    bool output_flag;
    bool log_to_console;
    HighsInt log_dev_level;
  };
  std::unique_ptr<LogData> log_data;
  HighsLogOptions log_options;

  bool use_original_HFactor_logic;
  bool debug_report_;
  HighsInt basis_matrix_limit_size;
  HighsInt update_method;

  // Internal timing
  HighsTimer* build_timer_;

  // Working buffer
  HighsInt nwork;
  vector<HighsInt> iwork;
  vector<double> dwork;

  // Basis matrix
  vector<HighsInt> b_var;  // Temp
  vector<HighsInt> b_start;
  vector<HighsInt> b_index;
  vector<double> b_value;

  // Permutation
  vector<HighsInt> permute;

  // Kernel matrix
  vector<HighsInt> mc_var;  // Temp
  vector<HighsInt> mc_start;
  vector<HighsInt> mc_count_a;
  vector<HighsInt> mc_count_n;
  vector<HighsInt> mc_space;
  vector<HighsInt> mc_index;
  vector<double> mc_value;
  vector<double> mc_min_pivot;

  // Row wise kernel matrix
  vector<HighsInt> mr_start;
  vector<HighsInt> mr_count;
  vector<HighsInt> mr_space;
  vector<HighsInt> mr_count_before;
  vector<HighsInt> mr_index;

  // Kernel column buffer
  vector<HighsInt> mwz_column_index;
  vector<HighsBool> mwz_column_mark;
  vector<double> mwz_column_array;

  // Count link list
  vector<HighsInt> col_link_first;
  vector<HighsInt> col_link_next;
  vector<HighsInt> col_link_last;

  vector<HighsInt> row_link_first;
  vector<HighsInt> row_link_next;
  vector<HighsInt> row_link_last;

  // Factor L
  vector<HighsInt> l_pivot_lookup;
  vector<HighsInt> l_pivot_index;

  vector<HighsInt> l_start;
  vector<HighsInt> l_index;
  vector<double> l_value;
  vector<HighsInt> lr_start;
  vector<HighsInt> lr_index;
  vector<double> lr_value;

  // Factor U
  vector<HighsInt> u_pivot_lookup;
  vector<HighsInt> u_pivot_index;
  vector<double> u_pivot_value;

  HighsInt u_merit_x;  // Only in PF and MPF
  HighsInt u_total_x;  // Only in PF and MPF
  vector<HighsInt> u_start;
  vector<HighsInt> u_last_p;
  vector<HighsInt> u_index;
  vector<double> u_value;
  vector<HighsInt> ur_start;
  vector<HighsInt> ur_lastp;
  vector<HighsInt> ur_space;
  vector<HighsInt> ur_index;
  vector<double> ur_value;

  // Update buffer
  vector<double> pf_pivot_value;
  vector<HighsInt> pf_pivot_index;
  vector<HighsInt> pf_start;
  vector<HighsInt> pf_index;
  vector<double> pf_value;

  HVector rhs_;

  // Implementation
  void buildSimple();
  //    void buildKernel();
  HighsInt buildKernel();
  void buildHandleRankDeficiency();
  void buildReportRankDeficiency();
  void buildMarkSingC();
  void buildFinish();
  void zeroCol(const HighsInt iCol, const double pivot_multiplier);
  void luClear();
  // Rebuild using refactor information
  HighsInt rebuild(HighsTimerClock* factor_timer_clock_pointer);

  // Action to take when pointers to the A matrix are no longer valid
  void invalidAMatrixAction();

  void reportIntVector(const std::string& name,
                       const vector<HighsInt>& entry) const;
  void reportDoubleVector(const std::string& name,
                          const vector<double>& entry) const;

  void ftranL(HVector& vector, const double expected_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void btranL(HVector& vector, const double expected_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void ftranU(HVector& vector, const double expected_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void btranU(HVector& vector, const double expected_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;

  void ftranFT(HVector& vector) const;
  void btranFT(HVector& vector) const;
  void ftranPF(HVector& vector) const;
  void btranPF(HVector& vector) const;
  void ftranMPF(HVector& vector) const;
  void btranMPF(HVector& vector) const;
  void ftranAPF(HVector& vector) const;
  void btranAPF(HVector& vector) const;

  void updateCFT(HVector* aq, HVector* ep, HighsInt* iRow);
  void updateFT(HVector* aq, HVector* ep, HighsInt iRow);
  void updatePF(HVector* aq, HighsInt iRow, HighsInt* hint);
  void updateMPF(HVector* aq, HVector* ep, HighsInt iRow, HighsInt* hint);
  void updateAPF(HVector* aq, HVector* ep, HighsInt iRow);

  /**
   * Local in-line functions
   */
  void colInsert(const HighsInt iCol, const HighsInt iRow, const double value) {
    const HighsInt iput = mc_start[iCol] + mc_count_a[iCol]++;
    mc_index[iput] = iRow;
    mc_value[iput] = value;
  }
  void colStoreN(const HighsInt iCol, const HighsInt iRow, const double value) {
    const HighsInt iput =
        mc_start[iCol] + mc_space[iCol] - (++mc_count_n[iCol]);
    mc_index[iput] = iRow;
    mc_value[iput] = value;
  }
  void colFixMax(const HighsInt iCol) {
    double max_value = 0;
    for (HighsInt k = mc_start[iCol]; k < mc_start[iCol] + mc_count_a[iCol];
         k++)
      max_value = max(max_value, fabs(mc_value[k]));
    // Record the minimum value for an acceptable pivot in this
    // column. Note that it can't just be max_value * pivot_threshold,
    // since this may be less than pivot_tolerance
    // (kDefaultPivotTolerance = 1e-10), possibly allowing a pivot
    // less than pivot_tolerance to be chosen rather than larger
    // values that satisfy the tolerance.
    //
    // Exposed by #3352 when the two entries in the pivotal column
    // were -8.1e-11 and -1.6e-10. With mc_min_pivot[iCol] = 1.6e-11,
    // -8.1e-11 was chosen as a pivot. However, it satisfies 8.1e-11 =
    // abs_pivot_value < pivot_tolerance so the singular pivot is
    // deferred.
    //
    // When singularity is detected, zeroCol is used to zero the
    // pivotal column, on the assumption that all entries are less
    // than the pivot (otherwise a larger pivot should have been
    // chosen). This triggers the original assert(abs_value <
    // pivot_tolerance), when abs_value = 1.6e-10. This is actually
    // conservative, since singularity should only be triggered if
    // there is no entry in the pivotal column >= pivot_tolerance, so
    // every pivotal column entry should be less than or equal to
    // abs_pivot_value.
    mc_min_pivot[iCol] = max(max_value * pivot_threshold, pivot_tolerance);
  }

  double colDelete(const HighsInt iCol, const HighsInt iRow) {
    HighsInt idel = mc_start[iCol];
    HighsInt imov = idel + (--mc_count_a[iCol]);
    while (mc_index[idel] != iRow) idel++;
    double pivot_multiplier = mc_value[idel];
    mc_index[idel] = mc_index[imov];
    mc_value[idel] = mc_value[imov];
    return pivot_multiplier;
  }

  void rowInsert(const HighsInt iCol, const HighsInt iRow) {
    HighsInt iput = mr_start[iRow] + mr_count[iRow]++;
    mr_index[iput] = iCol;
  }

  void rowDelete(const HighsInt iCol, const HighsInt iRow) {
    HighsInt idel = mr_start[iRow];
    HighsInt imov = idel + (--mr_count[iRow]);
    while (mr_index[idel] != iCol) idel++;
    mr_index[idel] = mr_index[imov];
  }

  void clinkAdd(const HighsInt index, const HighsInt count) {
    const HighsInt mover = col_link_first[count];
    col_link_last[index] = -2 - count;
    col_link_next[index] = mover;
    col_link_first[count] = index;
    if (mover >= 0) col_link_last[mover] = index;
  }

  void clinkDel(const HighsInt index) {
    const HighsInt xlast = col_link_last[index];
    const HighsInt xnext = col_link_next[index];
    if (xlast >= 0)
      col_link_next[xlast] = xnext;
    else
      col_link_first[-xlast - 2] = xnext;
    if (xnext >= 0) col_link_last[xnext] = xlast;
  }

  void rlinkAdd(const HighsInt index, const HighsInt count) {
    const HighsInt mover = row_link_first[count];
    row_link_last[index] = -2 - count;
    row_link_next[index] = mover;
    row_link_first[count] = index;
    if (mover >= 0) row_link_last[mover] = index;
  }

  void rlinkDel(const HighsInt index) {
    const HighsInt xlast = row_link_last[index];
    const HighsInt xnext = row_link_next[index];
    if (xlast >= 0)
      row_link_next[xlast] = xnext;
    else
      row_link_first[-xlast - 2] = xnext;
    if (xnext >= 0) row_link_last[xnext] = xlast;
  }
  friend class HSimplexNla;
};

#endif /* HFACTOR_H_ */
