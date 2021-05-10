/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDualRow.h
 * @brief Dual simplex ratio test for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HDUALROW_H_
#define SIMPLEX_HDUALROW_H_

#include <set>
#include <vector>

#include "lp_data/HighsModelObject.h"

class HVector;
const double initial_total_change = 1e-12;
const double initial_remain_theta = 1e100;
const double max_select_theta = 1e18;

/**
 * @brief Dual simplex ratio test for HiGHS
 *
 * Performs the dual bound-flipping ratio test and some update
 * dual/flip tasks
 */
class HDualRow {
 public:
  HDualRow(HighsModelObject& hmo) : workHMO(hmo) {}

  /**
   * @brief Calls setupSlice to set up the packed indices and values for
   * the dual ratio test
   */
  void setup();

  /**
   * @brief Set up the packed indices and values for the dual ratio test
   *
   * Done either for the whole pivotal row (see HDualRow::setup), or
   * just for a slice (see HDual::initSlice)
   */
  void setupSlice(int size  //!< Dimension of slice
  );
  /**
   * @brief Clear the packed data by zeroing packCount and workCount
   */
  void clear();

  /**
   * @brief Pack the indices and values for the row.
   *
   * Offset of numCol is used when packing row_ep
   */
  void chooseMakepack(const HVector* row,  //!< Row to be packed
                      const int offset     //!< Offset for indices
  );
  /**
   * @brief Determine the possible variables - candidates for CHUZC
   *
   * TODO: Check with Qi what this is doing
   */
  void choosePossible();

  /**
   * @brief Join pack of possible candidates in this row with possible
   * candidates in otherRow
   */
  void chooseJoinpack(const HDualRow* otherRow  //!< Other row to join with this
  );
  /**
   * @brief Chooses the entering variable via BFRT and EXPAND
   *
   * Can fail when there are excessive dual values due to EXPAND
   * perturbation not being relatively too small
   */
  bool chooseFinal();

  /**
   * @brief Identifies the groups of degenerate nodes in BFRT after a
   * heap sort of ratios
   */
  bool chooseFinalWorkGroupQuad();
  bool chooseFinalWorkGroupHeap();

  void chooseFinalLargeAlpha(
      int& breakIndex, int& breakGroup, int pass_workCount,
      const std::vector<std::pair<int, double>>& pass_workData,
      const std::vector<int>& pass_workGroup);

  void reportWorkDataAndGroup(
      const std::string message, const int reportWorkCount,
      const std::vector<std::pair<int, double>>& reportWorkData,
      const std::vector<int>& reportWorkGroup);
  bool compareWorkDataAndGroup();

  /**
   * @brief Update bounds when flips have occurred, and accumulate the
   * RHS for the FTRAN required to update the primal values after BFRT
   */
  void updateFlip(HVector* bfrtColumn  //!< RHS for FTRAN BFRT
  );
  /**
   * @brief Update the dual values
   */
  void updateDual(
      double theta  //!< Multiple of pivotal row to add int to duals
                    //      int columnOut  //!< Index of leaving column
  );
  /**
   * @brief Create a list of nonbasic free columns
   */
  void createFreelist();

  /**
   * @brief Set a value of nonbasicMove for all free columns to
   * prevent their dual values from being changed
   */
  void createFreemove(HVector* row_ep  //!< Row of \f$B^{-1}\f$ to be used to
                                       //!< compute pivotal row entry
  );
  /**
   * @brief Reset the nonbasicMove values for free columns
   */
  void deleteFreemove();

  /**
   * @brief Delete the list of nonbasic free columns
   */
  void deleteFreelist(int iColumn  //!< Index of column to remove from Freelist
  );

  /**
   * @brief Compute (contribution to) the Devex weight
   */
  void computeDevexWeight(const int slice = -1);

  HighsModelObject& workHMO;  //!< Local copy of pointer to model
  int workSize = -1;  //!< Size of the HDualRow slice: Initialise it here to
                      //!< avoid compiler warning
  const int* workNumTotPermutation;  //!< Pointer to model->numTotPermutation();
  const int* workMove;     //!< Pointer to workHMO.simplex_basis_.nonbasicMove_;
  const double* workDual;  //!< Pointer to workHMO.simplex_info_.workDual_;
  const double* workRange;  //!< Pointer to workHMO.simplex_info_.workRange_;
  const int*
      work_devex_index;  //!< Pointer to workHMO.simplex_info_.devex_index;

  // Freelist:
  std::set<int> freeList;  //!< Freelist itself

  // packed data:
  int packCount;                  //!< number of packed indices/values
  std::vector<int> packIndex;     //!< Packed indices
  std::vector<double> packValue;  //!< Packed values

  // (Local) value of computed weight
  double computed_edge_weight;

  double workDelta;  //!< Local copy of dual.deltaPrimal
  double workAlpha;  //!< Original copy of pivotal computed row-wise
  double workTheta;  //!< Original copy of dual step workDual[workPivot] /
                     //!< workAlpha;
  int workPivot;     //!< Index of the column entering the basis
  int workCount;     //!< Number of BFRT flips

  std::vector<std::pair<int, double>>
      workData;  //!< Index-Value pairs for ratio test
  std::vector<int>
      workGroup;  //!< Pointers into workData for degenerate nodes in BFRT

  // Independent identifiers for heap-based sort in BFRT
  int alt_workCount;
  std::vector<std::pair<int, double>> original_workData;
  std::vector<std::pair<int, double>> sorted_workData;
  std::vector<int> alt_workGroup;

  HighsSimplexAnalysis* analysis;
};

#endif /* SIMPLEX_HDUALROW_H_ */
