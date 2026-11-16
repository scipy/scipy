using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;

// mcs -out:highscslib.dll -t:library highs_csharp_api.cs -unsafe

namespace Highs {
public enum HighsStatus
{
    kError = -1,
    kOk,
    kWarning
}

public enum HighsMatrixFormat
{
    kColwise = 1,
    kRowwise
}

public enum HessianFormat
{
    kTriangular = 1,
    kSquare
}

public enum HighsBasisStatus
{
    kLower = 0,
    kBasic,
    kUpper,
    kZero,
    kNonbasic
}

public enum HighsObjectiveSense
{
    kMinimize = 1,
    kMaximize = -1
}

public enum HighsModelStatus
{
    kNotset = 0,
    kLoadError,
    kModelError,
    kPresolveError,
    kSolveError,
    kPostsolveError,
    kModelEmpty,
    kOptimal,
    kInfeasible,
    kUnboundedOrInfeasible,
    kUnbounded,
    kObjectiveBound,
    kObjectiveTarget,
    kTimeLimit,
    kIterationLimit,
    kUnknown,
    kSolutionLimit,
    kInterrupt,
    kMemoryLimit,
    kHighsInterrupt
}

public enum HighsIntegrality
{
    kContinuous = 0,
    kInteger = 1,
    kSemiContinuous = 2,
    kSemiInteger = 3,
    kImplicitInteger = 4,
}

public class HighsModel
{
    public HighsObjectiveSense sense;
    public double[] colcost;
    public double offset;
    public double[] collower;
    public double[] colupper;
    public double[] rowlower;
    public double[] rowupper;
    public HighsMatrixFormat a_format;
    public int[] astart;
    public int[] aindex;
    public double[] avalue;
    public int[] highs_integrality;

    public HighsModel()
    {

    }

    public HighsModel(double[] colcost, double[] collower, double[] colupper, double[] rowlower, double[] rowupper,
    int[] astart, int[] aindex, double[] avalue, int[] highs_integrality = null, double offset = 0, HighsMatrixFormat a_format = HighsMatrixFormat.kColwise, HighsObjectiveSense sense = HighsObjectiveSense.kMinimize)
    {
        this.colcost = colcost;
        this.collower = collower;
        this.colupper = colupper;
        this.rowlower = rowlower;
        this.rowupper = rowupper;
        this.astart = astart;
        this.aindex = aindex;
        this.avalue = avalue;
        this.offset = offset;
        this.a_format = a_format;
        this.sense = sense;
        this.highs_integrality = highs_integrality;
    }
}

public class HighsHessian
{
    public HessianFormat q_format;
    public int dim;
    public int[] qstart;
    public int[] qindex;
    public double[] qvalue;

    public HighsHessian()
    {

    }

    public HighsHessian(int dim, int[] qstart, int[] qindex, double[] qvalue, HessianFormat q_format = HessianFormat.kTriangular)
    {
        this.dim = dim;
        this.qstart = qstart;
        this.qindex = qindex;
        this.qvalue = qvalue;
        this.q_format = q_format;
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

public class HighsLpSolver : IDisposable
{
    private IntPtr highs;

    private bool _disposed;

    private const string highslibname = "highs";

    [DllImport(highslibname)]
    private static extern int Highs_call(
        Int32 numcol,
        Int32 numrow,
        Int32 numnz,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        double[] colvalue,
        double[] coldual,
        double[] rowvalue,
        double[] rowdual,
        int[] colbasisstatus,
        int[] rowbasisstatus,
        ref int modelstatus);

    [DllImport(highslibname)]
    private static extern IntPtr Highs_create();

    [DllImport(highslibname)]
    private static extern void Highs_destroy(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_run(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_readModel(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_writeModel(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_writePresolvedModel(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_writeSolutionPretty(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern double Highs_getInfinity(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_passLp(
        IntPtr highs,
        int numcol,
        int numrow,
        int numnz,
        int aformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue);

    [DllImport(highslibname)]
    private static extern int Highs_passMip(
        IntPtr highs,
        int numcol,
        int numrow,
        int numnz,
        int aformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        int[] highs_integrality);

    [DllImport(highslibname)]
    private static extern int Highs_passModel(
        IntPtr highs,
        int numcol,
        int numrow,
        int numnz,
        int qnumnz,
        int aformat,
        int qformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        int[] qstart,
        int[] qindex,
        double[] qvalue,
        int[] highs_integrality);

    [DllImport(highslibname)]
    private static extern int Highs_passHessian(
        IntPtr highs,
        int dim,
        int numnz,
        int q_format,
        int[] qstart,
        int[] qindex,
        double[] qvalue);

    [DllImport(highslibname)]
    private static extern int Highs_setOptionValue(IntPtr highs, string option, string value);

    [DllImport(highslibname)]
    private static extern int Highs_setBoolOptionValue(IntPtr highs, string option, int value);

    [DllImport(highslibname)]
    private static extern int Highs_setIntOptionValue(IntPtr highs, string option, int value);

    [DllImport(highslibname)]
    private static extern int Highs_setDoubleOptionValue(IntPtr highs, string option, double value);

    [DllImport(highslibname)]
    private static extern int Highs_setStringOptionValue(IntPtr highs, string option, string value);

    [DllImport(highslibname)]
    private static extern int Highs_getBoolOptionValue(IntPtr highs, string option, out int value);

    [DllImport(highslibname)]
    private static extern int Highs_getIntOptionValue(IntPtr highs, string option, out int value);

    [DllImport(highslibname)]
    private static extern int Highs_getDoubleOptionValue(IntPtr highs, string option, out double value);

    [DllImport(highslibname)]
    private static extern int Highs_getStringOptionValue(IntPtr highs, string option, [Out] StringBuilder value);

    [DllImport(highslibname)]
    private static extern int Highs_getSolution(IntPtr highs, double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual);

    [DllImport(highslibname)]
    private static extern int Highs_getNumCol(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getNumRow(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getNumNz(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getHessianNumNz(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getBasis(IntPtr highs, int[] colstatus, int[] rowstatus);    

    [DllImport(highslibname)]
    private static extern int Highs_setBasis(IntPtr highs, int[] col_status, int[] row_status);

    [DllImport(highslibname)]
    private static extern int Highs_setLogicalBasis(IntPtr highs);

    [DllImport(highslibname)]
    private static extern double Highs_getObjectiveValue(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getIterationCount(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getModelStatus(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_addRow(IntPtr highs, double lower, double upper, int num_new_nz, int[] indices, double[] values);

    [DllImport(highslibname)]
    private static extern int Highs_addRows(
        IntPtr highs,
        int num_new_row,
        double[] lower,
        double[] upper,
        int num_new_nz,
        int[] starts,
        int[] indices,
        double[] values);

    [DllImport(highslibname)]
    private static extern int Highs_addCol(
        IntPtr highs,
        double cost,
        double lower,
        double upper,
        int num_new_nz,
        int[] indices,
        double[] values);

    [DllImport(highslibname)]
    private static extern int Highs_addCols(
        IntPtr highs,
        int num_new_col,
        double[] costs,
        double[] lower,
        double[] upper,
        int num_new_nz,
        int[] starts,
        int[] indices,
        double[] values);

    [DllImport(highslibname)]
    private static extern int Highs_changeObjectiveSense(IntPtr highs, int sense);

    [DllImport(highslibname)]
    private static extern int Highs_changeColCost(IntPtr highs, int col, double cost);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsCostBySet(IntPtr highs, int num_set_entries, int[] set, double[] cost);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsCostByMask(IntPtr highs, int[] mask, double[] cost);

    [DllImport(highslibname)]
    private static extern int Highs_changeColBounds(IntPtr highs, int col, double lower, double upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsBoundsByRange(IntPtr highs, int from_col, int to_col, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsBoundsBySet(IntPtr highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsBoundsByMask(IntPtr highs, int[] mask, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeRowBounds(IntPtr highs, int row, double lower, double upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeRowsBoundsByRange(IntPtr highs, int from_row, int to_row, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeRowsBoundsBySet(IntPtr highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeRowsBoundsByMask(IntPtr highs, int[] mask, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsIntegralityByRange(IntPtr highs, int from_col, int to_col, int[] integrality);

    [DllImport(highslibname)]
    private static extern int Highs_changeCoeff(IntPtr highs, int row, int col, double value);

    [DllImport(highslibname)]
    private static extern int Highs_deleteColsByRange(IntPtr highs, int from_col, int to_col);

    [DllImport(highslibname)]
    private static extern int Highs_deleteColsBySet(IntPtr highs, int num_set_entries, int[] set);

    [DllImport(highslibname)]
    private static extern int Highs_deleteColsByMask(IntPtr highs, int[] mask);

    [DllImport(highslibname)]
    private static extern int Highs_deleteRowsByRange(IntPtr highs, int from_row, int to_row);

    [DllImport(highslibname)]
    private static extern int Highs_deleteRowsBySet(IntPtr highs, int num_set_entries, int[] set);

    [DllImport(highslibname)]
    private static extern int Highs_deleteRowsByMask(IntPtr highs, int[] mask);

    [DllImport(highslibname)]
    private static extern int Highs_getDoubleInfoValue(IntPtr highs, string info, out double value);

    [DllImport(highslibname)]
    private static extern int Highs_getIntInfoValue(IntPtr highs, string info, out int value);

    [DllImport(highslibname)]
    private static extern int Highs_getInt64InfoValue(IntPtr highs, string info, out long value);

    [DllImport(highslibname)]
    private static extern int Highs_changeObjectiveOffset(IntPtr highs, double offset);

    [DllImport(highslibname)]
    private static extern int Highs_addLinearObjective(IntPtr highs, double weight, double offset, double[] coefficients, double abs_tolerance, double rel_tolerance, int priority);

    [DllImport(highslibname)]
    private static extern int Highs_clearLinearObjectives(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_setSolution(IntPtr highs, double[] col_value, double[] row_value, double[] col_dual, double[] row_dual);

    [DllImport(highslibname)]
    private static extern int Highs_setSparseSolution(IntPtr highs, int num_entries, int[] index, double[] value);

    [DllImport(highslibname)]
    private static extern int Highs_getColsByRange(
        IntPtr highs,
        int from_col,
        int to_col,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getColsBySet(
        IntPtr highs,
        int num_set_entries,
        int[] set,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getColsByMask(
        IntPtr highs,
        int[] mask,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getRowsByRange(
        IntPtr highs,
        int from_row,
        int to_row,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getRowsBySet(
        IntPtr highs,
        int num_set_entries,
        int[] set,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getRowsByMask(
        IntPtr highs,
        int[] mask,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getBasicVariables(IntPtr highs, int[] basic_variables);

    [DllImport(highslibname)]
    private static extern int Highs_getBasisInverseRow(IntPtr highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getBasisInverseCol(IntPtr highs, int col, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getBasisSolve(
        IntPtr highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_num_nz,
        int[] solution_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getBasisTransposeSolve(
        IntPtr highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_nz,
        int[] solution_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getReducedRow(IntPtr highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getReducedColumn(IntPtr highs, int col, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [DllImport(highslibname)]
    private static extern int Highs_clearModel(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_clearSolver(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_passColName(IntPtr highs, int col, string name);

    [DllImport(highslibname)]
    private static extern int Highs_passRowName(IntPtr highs, int row, string name);

    [DllImport(highslibname)]
    private static extern int Highs_writeOptions(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_writeOptionsDeviations(IntPtr highs, string filename);

    public static HighsStatus call(HighsModel model, ref HighsSolution sol, ref HighsBasis bas, ref HighsModelStatus modelstatus)
    {
        int nc = model.colcost.Length;
        int nr = model.rowlower.Length;
        int nnz = model.avalue.Length;

        int[] colbasstat = new int[nc];
        int[] rowbasstat = new int[nr];

        int modelstate = 0;

        HighsStatus status = (HighsStatus)HighsLpSolver.Highs_call(
            nc,
            nr,
            nnz,
            model.colcost,
            model.collower,
            model.colupper,
            model.rowlower,
            model.rowupper,
            model.astart,
            model.aindex,
            model.avalue,
            sol.colvalue,
            sol.coldual,
            sol.rowvalue,
            sol.rowdual,
            colbasstat,
            rowbasstat,
            ref modelstate);

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
        this.Dispose(false);
    }

    public void Dispose()
    {
        this.Dispose(true);
        GC.SuppressFinalize(this);
    }

    protected virtual void Dispose(bool disposing)
    {
        if (this._disposed)
        {
            return;
        }

        HighsLpSolver.Highs_destroy(this.highs);
        this._disposed = true;
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

    public HighsStatus writePresolvedModel(string filename)
    {
        return (HighsStatus)HighsLpSolver.Highs_writePresolvedModel(this.highs, filename);
    }

    public HighsStatus writeSolutionPretty(string filename)
    {
        return (HighsStatus)HighsLpSolver.Highs_writeSolutionPretty(this.highs, filename);
    }

    public double getInfinity()
    {
        return HighsLpSolver.Highs_getInfinity(this.highs);
    }

    public HighsStatus passLp(HighsModel model)
    {
        return (HighsStatus)HighsLpSolver.Highs_passLp(
            this.highs,
            model.colcost.Length,
            model.rowlower.Length,
            model.avalue.Length,
            (int)model.a_format,
            (int)model.sense,
            model.offset,
            model.colcost,
            model.collower,
            model.colupper,
            model.rowlower,
            model.rowupper,
            model.astart,
            model.aindex,
            model.avalue);
    }

    public HighsStatus passMip(HighsModel model)
    {
        return (HighsStatus)HighsLpSolver.Highs_passMip(
            this.highs,
            model.colcost.Length,
            model.rowlower.Length,
            model.avalue.Length,
            (int)model.a_format,
            (int)model.sense,
            model.offset,
            model.colcost,
            model.collower,
            model.colupper,
            model.rowlower,
            model.rowupper,
            model.astart,
            model.aindex,
            model.avalue,
            model.highs_integrality);
    }

    public HighsStatus passHessian(HighsHessian hessian)
    {
        return (HighsStatus)HighsLpSolver.Highs_passHessian(
            this.highs,
            hessian.dim,
            hessian.qvalue.Length,
            (int)hessian.q_format,
            hessian.qstart,
            hessian.qindex,
            hessian.qvalue);
    }

    public HighsStatus setOptionValue(string option, string value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setOptionValue(this.highs, option, value);
    }

    public HighsStatus setStringOptionValue(string option, string value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setStringOptionValue(this.highs, option, value);
    }

    public HighsStatus setBoolOptionValue(string option, int value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setBoolOptionValue(this.highs, option, value);
    }

    public HighsStatus setDoubleOptionValue(string option, double value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setDoubleOptionValue(this.highs, option, value);
    }

    public HighsStatus setIntOptionValue(string option, int value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setIntOptionValue(this.highs, option, value);
    }

    public HighsStatus getStringOptionValue(string option, out string value)
    {
        var stringBuilder = new StringBuilder();
        var result = (HighsStatus)HighsLpSolver.Highs_getStringOptionValue(this.highs, option, stringBuilder);
        value = stringBuilder.ToString();
        return result;
    }

    public HighsStatus getBoolOptionValue(string option, out int value)
    {
        return (HighsStatus)HighsLpSolver.Highs_getBoolOptionValue(this.highs, option, out value);
    }

    public HighsStatus getDoubleOptionValue(string option, out double value)
    {
        return (HighsStatus)HighsLpSolver.Highs_getDoubleOptionValue(this.highs, option, out value);
    }

    public HighsStatus getIntOptionValue(string option, out int value)
    {
        return (HighsStatus)HighsLpSolver.Highs_getIntOptionValue(this.highs, option, out value);
    }

    public int getNumCol()
    {
        return HighsLpSolver.Highs_getNumCol(this.highs);
    }

    public int getNumRow()
    {
        return HighsLpSolver.Highs_getNumRow(this.highs);
    }

    public int getNumNz()
    {
        return HighsLpSolver.Highs_getNumNz(this.highs);
    }

    public HighsSolution getSolution()
    {
        int nc = this.getNumCol();
        int nr = this.getNumRow();

        HighsSolution sol = new HighsSolution(nc, nr);
        HighsLpSolver.Highs_getSolution(this.highs, sol.colvalue, sol.coldual, sol.rowvalue, sol.rowdual);

        return sol;
    }

    public HighsBasis getBasis()
    {
        int nc = this.getNumCol();
        int nr = this.getNumRow();

        int[] colbasstat = new int[nc];
        int[] rowbasstat = new int[nr];

        HighsLpSolver.Highs_getBasis(this.highs, colbasstat, rowbasstat);
        HighsBasis bas = new HighsBasis(
            colbasstat.Select(x => (HighsBasisStatus)x).ToArray(),
            rowbasstat.Select(x => (HighsBasisStatus)x).ToArray());

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
        return (HighsStatus)HighsLpSolver.Highs_addCols(
            this.highs,
            costs.Length,
            costs,
            lower,
            upper,
            indices.Length,
            starts,
            indices,
            values);
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

    public HighsStatus changeRowsBoundsByRange(int from, int to, double[] lower, double[] upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeRowsBoundsByRange(this.highs, from, to, lower, upper);
    }

    public HighsStatus changeRowsBoundsBySet(int[] rows, double[] lower, double[] upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeRowsBoundsBySet(this.highs, rows.Length, rows, lower, upper);
    }

    public HighsStatus changeRowsBoundsByMask(bool[] mask, double[] lower, double[] upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeRowsBoundsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), lower, upper);
    }

    public HighsStatus changeColsIntegralityByRange(int from_col, int to_col, HighsIntegrality[] integrality)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeColsIntegralityByRange(this.highs, from_col, to_col, Array.ConvertAll(integrality, item => (int)item));
    }

    public HighsStatus changeCoeff(int row, int col, double value)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeCoeff(this.highs, row, col, value);
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

    delegate int HighsGetInfoDelegate<TValue>(IntPtr highs, string infoName, out TValue output);

    private TValue GetValueOrFallback<TValue>(HighsGetInfoDelegate<TValue> highsGetInfoDelegate, string infoName, TValue fallback)
    {
        try
        {
            var status = (HighsStatus)highsGetInfoDelegate(this.highs, infoName, out var value);
            if (status != HighsStatus.kOk)
            {
                return fallback;
            }

            return value;
        }
        catch
        {
            return fallback;
        }
    }

    /// <summary>
    /// Gets the current solution info.
    /// </summary>
    /// <returns>The <see cref="SolutionInfo"/>.</returns>
    public SolutionInfo getInfo()
    {
        // TODO: This object does not contian the "complete" info from the C api. Add further props, if you need them.
        var info = new SolutionInfo()
        {
            MipGap = this.GetValueOrFallback(HighsLpSolver.Highs_getDoubleInfoValue, "mip_gap", double.NaN),
            DualBound = this.GetValueOrFallback(HighsLpSolver.Highs_getDoubleInfoValue, "mip_dual_bound", double.NaN),
            ObjectiveValue = this.GetValueOrFallback(HighsLpSolver.Highs_getDoubleInfoValue, "objective_function_value", double.NaN),
            NodeCount = this.GetValueOrFallback(HighsLpSolver.Highs_getInt64InfoValue, "mip_node_count", 0L),
            IpmIterationCount = this.GetValueOrFallback(HighsLpSolver.Highs_getIntInfoValue, "ipm_iteration_count", 0),
            SimplexIterationCount = this.GetValueOrFallback(HighsLpSolver.Highs_getIntInfoValue, "simplex_iteration_count", 0),
            PdlpIterationCount = this.GetValueOrFallback(HighsLpSolver.Highs_getIntInfoValue, "pdlp_iteration_count", 0),
        };
        return info;
    }

    public HighsStatus changeObjectiveOffset(double offset)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeObjectiveOffset(this.highs, offset);
    }

    public HighsStatus addLinearObjective(double weight, double offset, double[] coefficients, double abs_tolerance, double rel_tolerance, int priority)
    {
        return (HighsStatus)HighsLpSolver.Highs_addLinearObjective(this.highs, weight, offset, coefficients, abs_tolerance, rel_tolerance, priority);
    }

    public HighsStatus clearLinearObjectives()
    {
        return (HighsStatus)HighsLpSolver.Highs_clearLinearObjectives(this.highs);
    }

    public HighsStatus setSolution(HighsSolution solution)
    {
        return (HighsStatus)HighsLpSolver.Highs_setSolution(this.highs, solution.colvalue, solution.rowvalue, solution.coldual, solution.rowdual);
    }

    public HighsStatus setBasis(HighsBasis basis)
    {
        int[] col_status = basis.colbasisstatus.Select(x => (int)x).ToArray();
        int[] row_status = basis.rowbasisstatus.Select(x => (int)x).ToArray();
        return (HighsStatus)HighsLpSolver.Highs_setBasis(this.highs, col_status, row_status);
    }

    public HighsStatus setLogicalBasis()
    {
        return (HighsStatus)HighsLpSolver.Highs_setLogicalBasis(this.highs);
    }

    /// <summary>Set a partial primal solution by passing values for a set of variables</summary>
    /// <param name="valuesByIndex">A dictionary that maps variable indices to variable values</param>
    /// <remarks>The sparse solution set by this function has values for a subset of the model's variables.
    /// For each entry in <paramref name="valuesByIndex"/>, the key identifies a variable by index, and
    /// the value indicates the variable's value in the sparse solution.</remarks>
    /// <returns>A <see cref="HighsStatus"/> constant indicating whether the call succeeded</returns>
    public HighsStatus setSparseSolution(IReadOnlyDictionary<int, double> valuesByIndex)
    {
        return (HighsStatus)Highs_setSparseSolution(this.highs, valuesByIndex.Count, valuesByIndex.Keys.ToArray(), valuesByIndex.Values.ToArray());
    }

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

    public HighsStatus clearModel()
    {
        return (HighsStatus)Highs_clearModel(this.highs);
    }

    public HighsStatus clearSolver()
    {
        return (HighsStatus)Highs_clearSolver(this.highs);
    }

    public HighsStatus passColName(int col, string name)
    {
        return (HighsStatus)Highs_passColName(this.highs, col, name);
    }

    public HighsStatus passRowName(int row, string name)
    {
        return (HighsStatus)Highs_passRowName(this.highs, row, name);
    }

    public HighsStatus writeOptions(string filename)
    {
        return (HighsStatus)Highs_writeOptions(this.highs, filename);
    }

    public HighsStatus writeOptionsDeviations(string filename)
    {
        return (HighsStatus)Highs_writeOptionsDeviations(this.highs, filename);
    }

    public HighsStatus getColsByRange(int from_col, int to_col, ref int num_col, double[] costs, double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value)
    {
        return (HighsStatus)Highs_getColsByRange(this.highs, from_col, to_col, ref num_col, costs, lower, upper, ref num_nz, matrix_start, matrix_index, matrix_value);
    }

    public HighsStatus getColsBySet(int num_set_entries, int[] set, ref int num_col, double[] costs, double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value)
    {
        return (HighsStatus)Highs_getColsBySet(this.highs, num_set_entries, set, ref num_col, costs, lower, upper, ref num_nz, matrix_start, matrix_index, matrix_value);
    }

    public HighsStatus getColsByMask(int[] mask, ref int num_col, double[] costs, double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value)
    {
        return (HighsStatus)Highs_getColsByMask(this.highs, mask, ref num_col, costs, lower, upper, ref num_nz, matrix_start, matrix_index, matrix_value);
    }

    public HighsStatus getRowsByRange(int from_row, int to_row, ref int num_row, double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value)
    {
        return (HighsStatus)Highs_getRowsByRange(this.highs, from_row, to_row, ref num_row, lower, upper, ref num_nz, matrix_start, matrix_index, matrix_value);
    }

    public HighsStatus getRowsBySet(int num_set_entries, int[] set, ref int num_row, double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value)
    {
        return (HighsStatus)Highs_getRowsBySet(this.highs, num_set_entries, set, ref num_row, lower, upper, ref num_nz, matrix_start, matrix_index, matrix_value);
    }

    public HighsStatus getRowsByMask(int[] mask, ref int num_row, double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value)
    {
        return (HighsStatus)Highs_getRowsByMask(this.highs, mask, ref num_row, lower, upper, ref num_nz, matrix_start, matrix_index, matrix_value);
    }    
}

/// <summary>
/// The solution info.
/// </summary>
public class SolutionInfo
{
    /// <summary>
    /// Gets or sets the simplex iteration count.
    /// </summary>
    public int SimplexIterationCount { get; set; }

    /// <summary>
    /// Gets or sets the Interior Point Method (IPM) iteration count.
    /// </summary>
    public int IpmIterationCount { get; set; }

    /// <summary>
    /// Gets or sets the PDLP iteration count.
    /// </summary>
    public int PdlpIterationCount { get; set; }

    /// <summary>
    /// Gets or sets the MIP gap.
    /// </summary>
    public double MipGap { get; set; }

    /// <summary>
    /// Gets or sets the best dual bound.
    /// </summary>
    public double DualBound { get; set; }

    /// <summary>
    /// Gets or sets the MIP node count.
    /// </summary>
    public long NodeCount { get; set; }

    /// <summary>
    /// Gets or sets the objective value.
    /// </summary>
    public double ObjectiveValue { get; set; }
}
}
