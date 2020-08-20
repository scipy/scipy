/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HFactor.h
 * @brief Basis matrix factorization, update and solves for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef HFACTOR_H_
#define HFACTOR_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsAnalysis.h"

using std::max;
using std::min;
using std::vector;

class HVector;

enum UPDATE_METHOD {
  UPDATE_METHOD_FT = 1,
  UPDATE_METHOD_PF = 2,
  UPDATE_METHOD_MPF = 3,
  UPDATE_METHOD_APF = 4
};
/**
 * Necessary threshholds for historical density to trigger
 * hyper-sparse TRANs,
 */
const double hyperFTRANL = 0.15;
const double hyperFTRANU = 0.10;
const double hyperBTRANL = 0.10;
const double hyperBTRANU = 0.15;
/**
 * Necessary threshhold for RHS density to trigger hyper-sparse TRANs,
 */
const double hyperCANCEL = 0.05;
/**
 * Threshhold for result density for it to be considered as
 * hyper-sparse - only for reporting
 */
const double hyperRESULT = 0.10;
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
  /**
   * @brief Copy problem size and pointers of constraint matrix, and set
   * up space for INVERT
   *
   * Copy problem size and pointers to coefficient matrix, allocate
   * working buffer for INVERT, allocate space for basis matrix, L, U
   * factor and Update buffer, allocated space for Markowitz matrices,
   * count-link-list, L factor and U factor
   */
  void setup(int numCol,            //!< Number of columns
             int numRow,            //!< Number of rows
             const int* Astart,     //!< Column starts of constraint matrix
             const int* Aindex,     //!< Row indices of constraint matrix
             const double* Avalue,  //!< Row values of constraint matrix
             int* baseIndex,        //!< Indices of basic variables
             int highs_debug_level = HIGHS_DEBUG_LEVEL_MIN,
             FILE* logfile = NULL, FILE* output = NULL,
             int message_level = ML_NONE,
             const bool use_original_HFactor_logic = true,
             int updateMethod =
                 UPDATE_METHOD_FT  //!< Default update method is Forrest Tomlin
  );

  /**
   * @brief Form \f$PBQ=LU\f$ for basis matrix \f$B\f$ or report degree of rank
   * deficiency.
   *
   * @return 0 if successful, otherwise rank_deficiency>0
   *
   */
  int build(HighsTimerClock* factor_timer_clock_pointer = NULL);

  /**
   * @brief Solve \f$B\mathbf{x}=\mathbf{b}\f$ (FTRAN)
   */
  void ftran(HVector& vector,            //!< RHS vector \f$\mathbf{b}\f$
             double historical_density,  //!< Historical density of the result
             HighsTimerClock* factor_timer_clock_pointer = NULL) const;

  /**
   * @brief Solve \f$B^T\mathbf{x}=\mathbf{b}\f$ (BTRAN)
   */
  void btran(HVector& vector,            //!< RHS vector \f$\mathbf{b}\f$
             double historical_density,  //!< Historical density of the result
             HighsTimerClock* factor_timer_clock_pointer = NULL) const;

  /**
   * @brief Update according to
   * \f$B'=B+(\mathbf{a}_q-B\mathbf{e}_p)\mathbf{e}_p^T\f$
   */
  void update(HVector* aq,  //!< Vector \f$B^{-1}\mathbf{a}_q\f$
              HVector* ep,  //!< Vector \f$B^{-T}\mathbf{e}_p\f$
              int* iRow,    //!< Index of pivotal row
              int* hint     //!< Reinversion status
  );

  /**
   * @brief Wall clock time for INVERT
   */
  double build_realTick;

  /**
   * @brief The synthetic clock for INVERT
   */
  double build_syntheticTick;

  // Rank deficiency information

  /**
   * @brief Degree of rank deficiency in \f$B\f$
   */
  int rank_deficiency;

  /**
   * @brief Rows not pivoted on
   */
  vector<int> noPvR;

  /**
   * @brief Columns not pivoted on
   */
  vector<int> noPvC;

  /**
   * @brief Gets baseIndex since it is private
   */
  const int* getBaseIndex() const { return baseIndex; }

  /**
   * @brief Gets Astart since it is private
   */
  const int* getAstart() const { return Astart; }

  /**
   * @brief Gets Aindex since it is private
   */
  const int* getAindex() const { return Aindex; }

  /**
   * @brief Gets Avalue since it is private
   */
  const double* getAvalue() const { return Avalue; }

  // Properties of data held in HFactor.h. To "have" them means that
  // they are assigned.
  int haveArrays;
  // The representation of B^{-1} corresponds to the current basis
  int haveInvert;
  // The representation of B^{-1} corresponds to the current basis and is fresh
  int haveFreshInvert;
  int basis_matrix_num_el = 0;
  int invert_num_el = 0;
  int kernel_dim = 0;
  int kernel_num_el = 0;

  /**
   * Data of the factor
   */

  // private:
  // Problem size, coefficient matrix and update method
  int numRow;
  int numCol;

 private:
  const int* Astart;
  const int* Aindex;
  const double* Avalue;
  int* baseIndex;
  int updateMethod;
  bool use_original_HFactor_logic;
  int highs_debug_level;
  FILE* logfile;
  FILE* output;
  int message_level;

  // Working buffer
  int nwork;
  vector<int> iwork;
  vector<double> dwork;

  // Basis matrix
  vector<int> Bstart;
  vector<int> Bindex;
  vector<double> Bvalue;

  // Permutation
  vector<int> permute;

  // Kernel matrix
  vector<int> MCstart;
  vector<int> MCcountA;
  vector<int> MCcountN;
  vector<int> MCspace;
  vector<int> MCindex;
  vector<double> MCvalue;
  vector<double> MCminpivot;

  // Row wise kernel matrix
  vector<int> MRstart;
  vector<int> MRcount;
  vector<int> MRspace;
  vector<int> MRcountb4;
  vector<int> MRindex;

  // Kernel column buffer
  vector<int> McolumnIndex;
  vector<char> McolumnMark;
  vector<double> McolumnArray;

  // Count link list
  vector<int> clinkFirst;
  vector<int> clinkNext;
  vector<int> clinkLast;

  vector<int> rlinkFirst;
  vector<int> rlinkNext;
  vector<int> rlinkLast;

  // Factor L
  vector<int> LpivotLookup;
  vector<int> LpivotIndex;

  vector<int> Lstart;
  vector<int> Lindex;
  vector<double> Lvalue;
  vector<int> LRstart;
  vector<int> LRindex;
  vector<double> LRvalue;

  // Factor U
  vector<int> UpivotLookup;
  vector<int> UpivotIndex;
  vector<double> UpivotValue;

  int UmeritX;
  int UtotalX;
  vector<int> Ustart;
  vector<int> Ulastp;
  vector<int> Uindex;
  vector<double> Uvalue;
  vector<int> URstart;
  vector<int> URlastp;
  vector<int> URspace;
  vector<int> URindex;
  vector<double> URvalue;

  // Update buffer
  vector<double> PFpivotValue;
  vector<int> PFpivotIndex;
  vector<int> PFstart;
  vector<int> PFindex;
  vector<double> PFvalue;

  // Implementation
  void buildSimple();
  //    void buildKernel();
  int buildKernel();
  void buildHandleRankDeficiency();
  void buildReportRankDeficiency();
  void buildMarkSingC();
  void buildFinish();

  void ftranL(HVector& vector, double historical_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void btranL(HVector& vector, double historical_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void ftranU(HVector& vector, double historical_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void btranU(HVector& vector, double historical_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;

  void ftranFT(HVector& vector) const;
  void btranFT(HVector& vector) const;
  void ftranPF(HVector& vector) const;
  void btranPF(HVector& vector) const;
  void ftranMPF(HVector& vector) const;
  void btranMPF(HVector& vector) const;
  void ftranAPF(HVector& vector) const;
  void btranAPF(HVector& vector) const;

  void updateCFT(HVector* aq, HVector* ep, int* iRow);
  void updateFT(HVector* aq, HVector* ep, int iRow);
  void updatePF(HVector* aq, int iRow, int* hint);
  void updateMPF(HVector* aq, HVector* ep, int iRow, int* hint);
  void updateAPF(HVector* aq, HVector* ep, int iRow);

  /**
   * Local in-line functions
   */
  void colInsert(const int iCol, const int iRow, const double value) {
    const int iput = MCstart[iCol] + MCcountA[iCol]++;
    MCindex[iput] = iRow;
    MCvalue[iput] = value;
  }
  void colStoreN(const int iCol, const int iRow, const double value) {
    const int iput = MCstart[iCol] + MCspace[iCol] - (++MCcountN[iCol]);
    MCindex[iput] = iRow;
    MCvalue[iput] = value;
  }
  void colFixMax(const int iCol) {
    double maxValue = 0;
    for (int k = MCstart[iCol]; k < MCstart[iCol] + MCcountA[iCol]; k++)
      maxValue = max(maxValue, fabs(MCvalue[k]));
    MCminpivot[iCol] = maxValue * 0.1;
  }

  double colDelete(const int iCol, const int iRow) {
    int idel = MCstart[iCol];
    int imov = idel + (--MCcountA[iCol]);
    while (MCindex[idel] != iRow) idel++;
    double pivotX = MCvalue[idel];
    MCindex[idel] = MCindex[imov];
    MCvalue[idel] = MCvalue[imov];
    return pivotX;
  }

  void rowInsert(const int iCol, const int iRow) {
    int iput = MRstart[iRow] + MRcount[iRow]++;
    MRindex[iput] = iCol;
  }

  void rowDelete(const int iCol, const int iRow) {
    int idel = MRstart[iRow];
    int imov = idel + (--MRcount[iRow]);
    while (MRindex[idel] != iCol) idel++;
    MRindex[idel] = MRindex[imov];
  }

  void clinkAdd(const int index, const int count) {
    const int mover = clinkFirst[count];
    clinkLast[index] = -2 - count;
    clinkNext[index] = mover;
    clinkFirst[count] = index;
    if (mover >= 0) clinkLast[mover] = index;
  }

  void clinkDel(const int index) {
    const int xlast = clinkLast[index];
    const int xnext = clinkNext[index];
    if (xlast >= 0)
      clinkNext[xlast] = xnext;
    else
      clinkFirst[-xlast - 2] = xnext;
    if (xnext >= 0) clinkLast[xnext] = xlast;
  }

  void rlinkAdd(const int index, const int count) {
    const int mover = rlinkFirst[count];
    rlinkLast[index] = -2 - count;
    rlinkNext[index] = mover;
    rlinkFirst[count] = index;
    if (mover >= 0) rlinkLast[mover] = index;
  }

  void rlinkDel(const int index) {
    const int xlast = rlinkLast[index];
    const int xnext = rlinkNext[index];
    if (xlast >= 0)
      rlinkNext[xlast] = xnext;
    else
      rlinkFirst[-xlast - 2] = xnext;
    if (xnext >= 0) rlinkLast[xnext] = xlast;
  }
};

#endif /* HFACTOR_H_ */
