using System;
using System.Linq;
using System.Runtime.InteropServices;

// mcs -out:highscslib.dll -t:library highs_csharp_api.cs -unsafe

public enum HighsStatus
{
   OK,
   Warning,
   Error
}

public enum HighsBasisStatus
{
   Lower,
   Basic,
   Upper,
   Zero,
   Nonbasic,
   Super
}

public enum HighsObjectiveSense
{
   Minimize = 1,
   Maximize = -1
}

public enum HighsModelStatus
{
   NOTSET,
   LOAD_ERROR,
   MODEL_ERROR,
   MODEL_EMPTY,
   PRESOLVE_ERROR,
   SOLVE_ERROR,
   POSTSOLVE_ERROR,
   PRIMAL_FEASIBLE,
   DUAL_FEASIBLE,
   PRIMAL_INFEASIBLE,
   PRIMAL_UNBOUNDED,
   OPTIMAL,
   REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND,
   REACHED_TIME_LIMIT,
   REACHED_ITERATION_LIMIT
}

public class HighsModel
{
   public double[] colcost;
   public double[] collower;
   public double[] colupper;
   public double[] rowlower;
   public double[] rowupper;
   public int[] astart;
   public int[] aindex;
   public double[] avalue;

   public HighsModel()
   {

   }

   public HighsModel(double[] colcost, double[] collower, double[] colupper, double[] rowlower, double[] rowupper,
   int[] astart, int[] aindex, double[] avalue)
   {
      this.colcost = colcost;
      this.collower = collower;
      this.colupper = colupper;
      this.rowlower = rowlower;
      this.rowupper = rowupper;
      this.astart = astart;
      this.aindex = aindex;
      this.avalue = avalue;
   }
}

public class HighsSolution
{
   public double[] colvalue;
   public double[] coldual;
   public double[] rowvalue;
   public double[] rowdual;

   public HighsSolution(int numcol, int numrow)
   {
      this.colvalue = new double[numcol];
      this.coldual = new double[numcol];
      this.rowvalue = new double[numrow];
      this.rowdual = new double[numrow];
   }

   public HighsSolution(double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual)
   {
      this.colvalue = colvalue;
      this.coldual = coldual;
      this.rowvalue = rowvalue;
      this.rowdual = rowdual;
   }
}

public class HighsBasis
{
   public HighsBasisStatus[] colbasisstatus;
   public HighsBasisStatus[] rowbasisstatus;

   public HighsBasis(int numcol, int numrow)
   {
      this.colbasisstatus = new HighsBasisStatus[numcol];
      this.rowbasisstatus = new HighsBasisStatus[numrow];
   }

   public HighsBasis(HighsBasisStatus[] colbasisstatus, HighsBasisStatus[] rowbasisstatus)
   {
      this.colbasisstatus = colbasisstatus;
      this.rowbasisstatus = rowbasisstatus;
   }
}

public unsafe class HighsLpSolver
{
   private void* highs;

   private const string highslibname = "highs.dll";

   [DllImport(highslibname)]
   private static extern int Highs_call(Int32 numcol, Int32 numrow, Int32 numnz, double[] colcost,
   double[] collower, double[] colupper, double[] rowlower, double[] rowupper, int[] astart, int[] aindex, double[] avalue,
   double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual, int[] colbasisstatus, int[] rowbasisstatus, ref int modelstatus);

   [DllImport(highslibname)]
   private static extern void* Highs_create();

   [DllImport(highslibname)]
   private static extern void Highs_destroy(void* highs);

   [DllImport(highslibname)]
   private static extern int Highs_run(void* highs);

   [DllImport(highslibname)]
   private static extern int Highs_readModel(void* highs, string filename);

   [DllImport(highslibname)]
   private static extern int Highs_writeModel(void* highs, string filename);

   [DllImport(highslibname)]
   private static extern int Highs_passLp(void* highs, int numcol, int numrow, int numnz, double[] colcost,
   double[] collower, double[] colupper, double[] rowlower, double[] rowupper, int[] astart, int[] aindex, double[] avalue);

   [DllImport(highslibname)]
   private static extern int Highs_setOptionValue(void* highs, string option, string value);

   [DllImport(highslibname)]
   private static extern void Highs_getSolution(void* highs, double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual);

   [DllImport(highslibname)]
   private static extern int Highs_getNumCols(void* highs);

   [DllImport(highslibname)]
   private static extern int Highs_getNumRows(void* highs);

   [DllImport(highslibname)]
   private static extern int Highs_getNumNz(void* highs);

   [DllImport(highslibname)]
   private static extern void Highs_getBasis(void* highs, int[] colstatus, int[] rowstatus);

   [DllImport(highslibname)]
   private static extern double Highs_getObjectiveValue(void* highs);

   [DllImport(highslibname)]
   private static extern int Highs_getIterationCount(void* highs);

   [DllImport(highslibname)]
   private static extern int Highs_getModelStatus(void* highs);

   [DllImport(highslibname)]
   private static extern int Highs_addRow(void* highs, double lower, double upper, int num_new_nz, int[] indices, double[] values);

   [DllImport(highslibname)]
   private static extern int Highs_addRows(void* highs, int num_new_row, double[] lower, double[] upper,
   int num_new_nz, int[] starts, int[] indices, double[] values);

   [DllImport(highslibname)]
   private static extern int Highs_addCol(void* highs, double cost, double lower, double upper,
   int num_new_nz, int[] indices, double[] values);

   [DllImport(highslibname)]
   private static extern int Highs_addCols(void* highs, int num_new_col, double[] costs, double[] lower, double[] upper,
   int num_new_nz, int[] starts, int[] indices, double[] values);

   [DllImport(highslibname)]
   private static extern int Highs_changeObjectiveSense(void* highs, int sense);

   [DllImport(highslibname)]
   private static extern int Highs_changeColCost(void* highs, int col, double cost);

   [DllImport(highslibname)]
   private static extern int Highs_changeColsCostBySet(void* highs, int num_set_entries, int[] set, double[] cost);

   [DllImport(highslibname)]
   private static extern int Highs_changeColsCostByMask(void* highs, int[] mask, double[] cost);

   [DllImport(highslibname)]
   private static extern int Highs_changeColBounds(void* highs, int col, double lower, double upper);

   [DllImport(highslibname)]
   private static extern int Highs_changeColsBoundsByRange(void* highs, int from_col, int to_col, double[] lower, double[] upper);

   [DllImport(highslibname)]
   private static extern int Highs_changeColsBoundsBySet(void* highs, int num_set_entries, int[] set, double[] lower, double[] upper);

   [DllImport(highslibname)]
   private static extern int Highs_changeColsBoundsByMask(void* highs, int[] mask, double[] lower, double[] upper);

   [DllImport(highslibname)]
   private static extern int Highs_changeRowBounds(void* highs, int row, double lower, double upper);

   [DllImport(highslibname)]
   private static extern int Highs_changeRowsBoundsBySet(void* highs, int num_set_entries, int[] set, double[] lower, double[] upper);

   [DllImport(highslibname)]
   private static extern int Highs_changeRowsBoundsByMask(void* highs, int[] mask, double[] lower, double[] upper);

   [DllImport(highslibname)]
   private static extern int Highs_deleteColsByRange(void* highs, int from_col, int to_col);

   [DllImport(highslibname)]
   private static extern int Highs_deleteColsBySet(void* highs, int num_set_entries, int[] set);

   [DllImport(highslibname)]
   private static extern int Highs_deleteColsByMask(void* highs, int[] mask);

   [DllImport(highslibname)]
   private static extern int Highs_deleteRowsByRange(void* highs, int from_row, int to_row);

   [DllImport(highslibname)]
   private static extern int Highs_deleteRowsBySet(void* highs, int num_set_entries, int[] set);

   [DllImport(highslibname)]
   private static extern int Highs_deleteRowsByMask(void* highs, int[] mask);

   [DllImport(highslibname)]
   private static extern int Highs_getColsByRange(void* highs, int from_col, int to_col, ref int num_col, double[] costs,
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport(highslibname)]
   private static extern int Highs_getColsBySet(void* highs, int num_set_entries, int[] set, ref int num_col, double[] costs,
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport(highslibname)]
   private static extern int Highs_getColsByMask(void* highs, int[] mask, ref int num_col, double[] costs,
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport(highslibname)]
   private static extern int Highs_getRowsByRange(void* highs, int from_row, int to_row, ref int num_row,
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport(highslibname)]
   private static extern int Highs_getRowsBySet(void* highs, int num_set_entries, int[] set, ref int num_row,
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport(highslibname)]
   private static extern int Highs_getRowsByMask(void* highs, int[] mask, ref int num_row,
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport(highslibname)]
   private static extern int Highs_getBasicVariables(void* highs, int[] basic_variables);

   [DllImport(highslibname)]
   private static extern int Highs_getBasisInverseRow(void* highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

   [DllImport(highslibname)]
   private static extern int Highs_getBasisInverseCol(void* highs, int col, double[] col_vector, ref int col_num_nz, int[] col_indices);

   [DllImport(highslibname)]
   private static extern int Highs_getBasisSolve(void* highs, double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices);

   [DllImport(highslibname)]
   private static extern int Highs_getBasisTransposeSolve(void* highs, double[] rhs, double[] solution_vector, ref int solution_nz, int[] solution_indices);

   [DllImport(highslibname)]
   private static extern int Highs_getReducedRow(void* highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

   [DllImport(highslibname)]
   private static extern int Highs_getReducedColumn(void* highs, int col, double[] col_vector, ref int col_num_nz, int[] col_indices);

   public static HighsStatus call(HighsModel model, ref HighsSolution sol, ref HighsBasis bas, ref HighsModelStatus modelstatus)
   {
      int nc = model.colcost.Length;
      int nr = model.rowlower.Length;
      int nnz = model.avalue.Length;

      int[] colbasstat = new int[nc];
      int[] rowbasstat = new int[nr];

      int modelstate = 0;

      HighsStatus status = (HighsStatus)HighsLpSolver.Highs_call(nc, nr, nnz, model.colcost, model.collower, model.colupper,
      model.rowlower, model.rowupper, model.astart, model.aindex, model.avalue,
      sol.colvalue, sol.coldual, sol.rowvalue, sol.rowdual, colbasstat, rowbasstat, ref modelstate);

      modelstatus = (HighsModelStatus)modelstate;

      bas.colbasisstatus = colbasstat.Select(x => (HighsBasisStatus)x).ToArray();
      bas.rowbasisstatus = rowbasstat.Select(x => (HighsBasisStatus)x).ToArray();

      return status;
   }

   public HighsLpSolver()
   {
      this.highs = HighsLpSolver.Highs_create();
   }

   ~HighsLpSolver()
   {
      HighsLpSolver.Highs_destroy(this.highs);
   }

   public HighsStatus run()
   {
      return (HighsStatus)HighsLpSolver.Highs_run(this.highs);
   }

   public HighsStatus readModel(string filename)
   {
      return (HighsStatus)HighsLpSolver.Highs_readModel(this.highs, filename);
   }

   public HighsStatus writeModel(string filename)
   {
      return (HighsStatus)HighsLpSolver.Highs_writeModel(this.highs, filename);
   }

   public HighsStatus passLp(HighsModel model)
   {
      return (HighsStatus)HighsLpSolver.Highs_passLp(this.highs, model.colcost.Length, model.rowlower.Length, model.avalue.Length,
      model.colcost, model.collower, model.colupper, model.rowlower, model.rowupper, model.astart, model.aindex, model.avalue);
   }

   public HighsStatus setOptionValue(string option, string value)
   {
      return (HighsStatus)HighsLpSolver.Highs_setOptionValue(this.highs, option, value);
   }

   public int getNumCols()
   {
      return HighsLpSolver.Highs_getNumCols(this.highs);
   }

   public int getNumRows()
   {
      return HighsLpSolver.Highs_getNumRows(this.highs);
   }

   public int getNumNz()
   {
      return HighsLpSolver.Highs_getNumNz(this.highs);
   }

   public HighsSolution getSolution()
   {
      int nc = this.getNumCols();
      int nr = this.getNumRows();

      HighsSolution sol = new HighsSolution(nc, nr);
      HighsLpSolver.Highs_getSolution(this.highs, sol.colvalue, sol.coldual, sol.rowvalue, sol.rowdual);
      return sol;
   }

   public HighsBasis getBasis()
   {
      int nc = this.getNumCols();
      int nr = this.getNumRows();

      int[] colbasstat = new int[nc];
      int[] rowbasstat = new int[nr];

      HighsLpSolver.Highs_getBasis(this.highs, colbasstat, rowbasstat);
      HighsBasis bas = new HighsBasis(colbasstat.Select(x => (HighsBasisStatus)x).ToArray(), rowbasstat.Select(x => (HighsBasisStatus)x).ToArray());

      return bas;
   }

   public double getObjectiveValue()
   {
      return HighsLpSolver.Highs_getObjectiveValue(this.highs);
   }

   public HighsModelStatus GetModelStatus()
   {
      return (HighsModelStatus)HighsLpSolver.Highs_getModelStatus(this.highs);
   }

   public int getIterationCount()
   {
      return HighsLpSolver.Highs_getIterationCount(this.highs);
   }

   public HighsStatus addRow(double lower, double upper, int[] indices, double[] values)
   {
      return (HighsStatus)HighsLpSolver.Highs_addRow(this.highs, lower, upper, indices.Length, indices, values);
   }

   public HighsStatus addRows(double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
   {
      return (HighsStatus)HighsLpSolver.Highs_addRows(this.highs, lower.Length, lower, upper, indices.Length, starts, indices, values);
   }

   public HighsStatus addCol(double cost, double lower, double upper, int[] indices, double[] values)
   {
      return (HighsStatus)HighsLpSolver.Highs_addCol(this.highs, cost, lower, upper, indices.Length, indices, values);
   }

   public HighsStatus addCols(double[] costs, double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
   {
      return (HighsStatus)HighsLpSolver.Highs_addCols(this.highs, costs.Length, costs, lower, upper, indices.Length, starts, indices, values);
   }

   public HighsStatus changeObjectiveSense(HighsObjectiveSense sense)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeObjectiveSense(this.highs, (int)sense);
   }

   public HighsStatus changeColCost(int col, double cost)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeColCost(this.highs, col, cost);
   }

   public HighsStatus changeColsCostBySet(int[] cols, double[] costs)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeColsCostBySet(this.highs, cols.Length, cols, costs);
   }

   public HighsStatus changeColsCostByMask(bool[] mask, double[] cost)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeColsCostByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), cost);
   }

   public HighsStatus changeColBounds(int col, double lower, double upper)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeColBounds(this.highs, col, lower, upper);
   }

   public HighsStatus changeColsBoundsByRange(int from, int to, double[] lower, double[] upper)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeColsBoundsByRange(this.highs, from, to, lower, upper);
   }

   public HighsStatus changeColsBoundsBySet(int[] cols, double[] lower, double[] upper)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeColsBoundsBySet(this.highs, cols.Length, cols, lower, upper);
   }

   public HighsStatus changeColsBoundsByMask(bool[] mask, double[] lower, double[] upper)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeColsBoundsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), lower, upper);
   }

   public HighsStatus changeRowBounds(int row, double lower, double upper)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeRowBounds(this.highs, row, lower, upper);
   }

   public HighsStatus changeRowsBoundsBySet(int[] rows, double[] lower, double[] upper)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeRowsBoundsBySet(this.highs, rows.Length, rows, lower, upper);
   }

   public HighsStatus changeRowsBoundsByMask(bool[] mask, double[] lower, double[] upper)
   {
      return (HighsStatus)HighsLpSolver.Highs_changeRowsBoundsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), lower, upper);
   }

   public HighsStatus deleteColsByRange(int from, int to)
   {
      return (HighsStatus)HighsLpSolver.Highs_deleteColsByRange(this.highs, from, to);
   }

   public HighsStatus deleteColsBySet(int[] cols)
   {
      return (HighsStatus)HighsLpSolver.Highs_deleteColsBySet(this.highs, cols.Length, cols);
   }

   public HighsStatus deleteColsByMask(bool[] mask)
   {
      return (HighsStatus)HighsLpSolver.Highs_deleteColsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray());
   }

   public HighsStatus deleteRowsByRange(int from, int to)
   {
      return (HighsStatus)HighsLpSolver.Highs_deleteRowsByRange(this.highs, from, to);
   }

   public HighsStatus deleteRowsBySet(int[] rows)
   {
      return (HighsStatus)HighsLpSolver.Highs_deleteRowsBySet(this.highs, rows.Length, rows);
   }

   public HighsStatus deleteRowsByMask(bool[] mask)
   {
      return (HighsStatus)HighsLpSolver.Highs_deleteRowsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray());
   }

   // int Highs_getColsByRange(void *highs, int from_col, int to_col, ref int num_col, double[] costs, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   // [DllImport(highslibname)]
   // int Highs_getColsBySet(void *highs, int num_set_entries, int[] set, ref int num_col, double[] costs, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   // [DllImport(highslibname)]
   // int Highs_getColsByMask(void *highs, int[] mask, ref int num_col, double[] costs, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   // [DllImport(highslibname)]
   // int Highs_getRowsByRange(void *highs, int from_row, int to_row, ref int num_row, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   // [DllImport(highslibname)]
   // int Highs_getRowsBySet(void *highs, int num_set_entries, int[] set, ref int num_row, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   // [DllImport(highslibname)]
   // int Highs_getRowsByMask(void *highs, int[] mask, ref int num_row, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   public HighsStatus getBasicVariables(ref int[] basic_variables)
   {
      return (HighsStatus)Highs_getBasicVariables(this.highs, basic_variables);
   }

   public HighsStatus getBasisInverseRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
   {
      return (HighsStatus)Highs_getBasisInverseRow(this.highs, row, row_vector, ref row_num_nz, row_indices);
   }

   public HighsStatus getBasisInverseCol(int col, double[] col_vector, ref int col_num_nz, int[] col_indices)
   {
      return (HighsStatus)Highs_getBasisInverseCol(this.highs, col, col_vector, ref col_num_nz, col_indices);
   }

   public HighsStatus getBasisSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
   {
      return (HighsStatus)Highs_getBasisSolve(this.highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
   }

   public HighsStatus getBasisTransposeSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
   {
      return (HighsStatus)Highs_getBasisTransposeSolve(this.highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
   }

   public HighsStatus getReducedRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
   {
      return (HighsStatus)Highs_getReducedRow(this.highs, row, row_vector, ref row_num_nz, row_indices);
   }

   public HighsStatus getReducedColumn(int col, double[] col_vector, ref int col_num_nz, int[] col_indices)
   {
      return (HighsStatus)Highs_getReducedColumn(this.highs, col, col_vector, ref col_num_nz, col_indices);
   }
}
