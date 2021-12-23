/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPresolve.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_PRESOLVE_H_
#define PRESOLVE_PRESOLVE_H_

#include <list>
#include <map>
#include <stack>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"
#include "presolve/HPreData.h"
#include "presolve/PresolveAnalysis.h"
#include "test/DevKkt.h"

using std::list;
using std::string;

enum class HighsPostsolveStatus {
  NotPresolved = -1,
  ReducedSolutionEmpty,
  ReducedSolutionDimenionsError,
  SolutionRecovered,
  LpOrPresolveObjectMissing,
  BasisError,
  NoPostsolve
};

enum class HighsPresolveStatus {
  NotPresolved = -1,
  NotReduced,
  Infeasible,
  Unbounded,
  Empty,
  Reduced,
  ReducedToEmpty,
  Timeout,
  NullError,
  OptionsError,
};

namespace presolve {

enum class Presolver {
  kMainEmpty,
  kMainRowSingletons,
  kMainForcing,
  kMainColSingletons,
  kMainDoubletonEq,
  kMainDominatedCols,
  kMainSingletonsOnly,
};

const std::map<Presolver, std::string> kPresolverNames{
    {Presolver::kMainEmpty, "Empty & fixed ()"},
    {Presolver::kMainRowSingletons, "Row singletons ()"},
    {Presolver::kMainForcing, "Forcing rows ()"},
    {Presolver::kMainColSingletons, "Col singletons ()"},
    {Presolver::kMainDoubletonEq, "Doubleton eq ()"},
    {Presolver::kMainDominatedCols, "Dominated Cols()"},
    {Presolver::kMainSingletonsOnly, "Singletons only()"}};

class Presolve : public HPreData {
 public:
  Presolve(HighsTimer& timer_ref) : timer(timer_ref) {}
  virtual ~Presolve() {}

  HighsPresolveStatus presolve();
  HighsPostsolveStatus postsolve(const HighsSolution& reduced_solution,
                                 const HighsBasis& reduced_basis,
                                 HighsSolution& recovered_solution,
                                 HighsBasis& recovered_basis);

  void setNumericalTolerances();
  void load(const HighsLp& lp);
  // todo: clear the public from below.
  string modelName;

  // Options
  std::vector<Presolver> order;

  int max_iterations = 0;

  void setTimeLimit(const double limit) {
    assert(limit < inf && limit > 0);
    timer.time_limit = limit;
  }

  int iPrint = 0;
  int message_level;
  FILE* output;

 private:
  int iKKTcheck = 0;
  int presolve(int print);

  const bool report_postsolve = false;

  double objShift;
  void initializeVectors();
  void setProblemStatus(const int s);
  void reportTimes();

  // new bounds on primal variables for implied free detection
  vector<double> implColLower;
  vector<double> implColUpper;
  vector<int> implColLowerRowIndex;
  vector<int> implColUpperRowIndex;

  vector<int> implRowDualLowerSingColRowIndex;
  vector<int> implRowDualUpperSingColRowIndex;

  // new bounds on row duals y_i
  vector<double> implRowDualLower;
  vector<double> implRowDualUpper;

  vector<double> implColDualLower;
  vector<double> implColDualUpper;
  vector<double> implRowValueLower;
  vector<double> implRowValueUpper;

  PresolveTimer timer;  // holds enum for main presolve rules

  enum stat {
    Unset = 0,
    Infeasible = 1,
    Unbounded = 2,
    Empty = 3,
    Optimal = 4,
    Reduced = 5,
    Timeout = 6,
  };

 private:
  bool hasChange = true;
  int status = 0;  // 0 is unassigned, see enum stat

  list<int> singRow;  // singleton rows
  list<int> singCol;  // singleton columns

  // original data
 public:
  vector<double> colCostOriginal;

 private:
  vector<double> rowLowerOriginal;
  vector<double> rowUpperOriginal;
  vector<double> colLowerOriginal;
  vector<double> colUpperOriginal;

  // functions
  void setPrimalValue(const int j, const double value);
  void checkForChanges(int iteration);
  void resizeProblem();
  void resizeImpliedBounds();

  // easy transformations
  void removeFixedCol(int j);
  void removeEmpty();
  void removeFixed();
  void removeEmptyRow(int i);
  void removeEmptyColumn(int j);
  void removeRow(int i);
  void checkBoundsAreConsistent();

  // singleton rows
  void removeRowSingletons();
  int getSingRowElementIndexInAR(int i);
  int getSingColElementIndexInA(int j);

  // forcing constraints
  void removeForcingConstraints();
  pair<double, double> getImpliedRowBounds(int row);
  void setVariablesToBoundForForcingRow(const int row, const bool isLower);
  void dominatedConstraintProcedure(const int i, const double g,
                                    const double h);

  // doubleton equations
  void removeDoubletonEquations();
  pair<int, int> getXYDoubletonEquations(const int row);
  void processRowDoubletonEquation(const int row, const int x, const int y,
                                   const double akx, const double aky,
                                   const double b);
  pair<double, double> getNewBoundsDoubletonConstraint(int row, int col, int j,
                                                       double aik, double aij);
  void UpdateMatrixCoeffDoubletonEquationXzero(const int i, const int x,
                                               const int y, const double aiy,
                                               const double akx,
                                               const double aky);
  void UpdateMatrixCoeffDoubletonEquationXnonZero(const int i, const int x,
                                                  const int y, const double aiy,
                                                  const double akx,
                                                  const double aky);

  // column singletons
  void removeColumnSingletons();
  bool removeIfImpliedFree(int col, int i, int k);
  void removeFreeColumnSingleton(const int col, const int row, const int k);
  void removeZeroCostColumnSingleton(const int col, const int row, const int k);
  bool removeColumnSingletonInDoubletonInequality(const int col, const int i,
                                                  const int k);
  void removeSecondColumnSingletonInDoubletonRow(const int j, const int i);
  pair<double, double> getBoundsImpliedFree(double lowInit, double uppInit,
                                            const int col, const int i,
                                            const int k);
  void removeImpliedFreeColumn(const int col, const int i, const int k);

  // dominated columns
  void removeDominatedColumns();
  void rowDualBoundsDominatedColumns();
  pair<double, double> getImpliedColumnBounds(int j);
  void removeIfWeaklyDominated(const int j, const double d, const double e);

  //    void findDuplicateRows();
  //    void findDuplicateColumns();
  //    void removeDuplicateRows(int i, int k, double v);
  //    int makeCheckForDuplicateRows(int k, int i, vector<double>& coeff,
  //    vector<int>& colIndex, double v, int whichIsFirst); void
  //    removeDuplicateColumns(int j,int k, double v); bool
  //    checkDuplicateRows(int i, int k) ;
  //	  bool checkDuplicateColumns(int i, int k) ;

  // old or test
  // void updateRemovedColRow(int dim);
  // void updateRowsByNZ();
  void testAnAR(int post);

  void countRemovedRows(PresolveRule rule);
  void countRemovedCols(PresolveRule rule);

  double tol = 0.0000001;
  const double default_primal_feasiblility_tolerance = 1e-7;
  const double default_dual_feasiblility_tolerance = 1e-7;
  const double default_small_matrix_value = 1e-9;
  double inconsistent_bounds_tolerance;
  double fixed_column_tolerance;
  double doubleton_equation_bound_tolerance;
  double doubleton_inequality_bound_tolerance;
  double presolve_small_matrix_value;
  double empty_row_bound_tolerance;
  double dominated_column_tolerance;
  double weakly_dominated_column_tolerance;

  // postsolve
  bool noPostSolve = false;

  void addChange(const PresolveRule type, const int row, const int col);
  void fillStackRowBounds(const int col);
  void setKKTcheckerData();

  void getBoundOnLByZj(const int row, const int j, double* lo, double* up,
                       const double colLow, const double colUpp);
  double getRowDualPost(const int row, const int col);
  double getColumnDualPost(const int col);
  string getDualsForcingRow(const int row, vector<int>& fRjs);
  void getDualsSingletonRow(const int row, const int col);
  void getDualsDoubletonEquation(const int row, const int col);
  void recordCounts(const string fileName);
  void trimA();

  void setBasisElement(const change c);

  // test basis matrix singularity
  //
  // public:
  //	vector<int> nbffull;
  //	vector<int> bindfull;
  //	void cmpNBF(int row, int col);
  //	void setNBFfullproblem(vector<int>& nbfFull, vector<int>& bnFull);
  //	int testBasisMatrixSingularity();
  //

  // Dev presolve
  // April 2020
  void reportDevMainLoop();
  void reportDevMidMainLoop();
  PresolveStats stats;
  int runPresolvers(const std::vector<Presolver>& order);

  void checkKkt(const bool final = false);
  dev_kkt_check::State initState(const bool intermediate = false);

  void caseTwoSingletonsDoubletonInequality(const int row, const int x,
                                            const int y);

  // August 2020
  void removeSingletonsOnly();
  bool isKnapsack(const int col) const;
  void removeKnapsack(const int col);
};

}  // namespace presolve

#endif /* PRESOLVE_HPRESOLVE_H_ */
