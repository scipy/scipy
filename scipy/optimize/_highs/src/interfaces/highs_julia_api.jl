struct HighsModel
   colcost
   collower
   colupper
   rowlower
   rowupper
   astart
   aindex
   avalue
end

struct HighsSolution
   colvalue
   coldual
   rowvalue
   rowdual
end

struct HighsBasis
   colbasisstatus
   rowbasisstatus
end

function Highs_call(model)
   n_col = convert(Cint, size(model.colcost, 1))
   n_row = convert(Cint, size(model.rowlower, 1))
   n_nz = convert(Cint, size(model.aindex, 1))

   colcost = convert(Array{Cdouble}, model.colcost)
   collower = convert(Array{Cdouble}, model.collower)
   colupper = convert(Array{Cdouble}, model.colupper)

   rowlower = convert(Array{Cdouble}, model.rowlower)
   rowupper = convert(Array{Cdouble}, model.rowupper)
   matstart = convert(Array{Cint}, model.astart)
   matindex = convert(Array{Cint}, model.aindex)
   matvalue = convert(Array{Cdouble}, model.avalue)

   solution = HighsSolution(Array{Cdouble, 1}(undef, n_col), Array{Cdouble, 1}(undef, n_col), Array{Cdouble, 1}(undef, n_row),  Array{Cdouble, 1}(undef, n_row))
   basis = HighsBasis(Array{Cint, 1}(undef, n_col), Array{Cint, 1}(undef, n_row)) 

   modelstatus = convert(Cint, 0)

   status = ccall((:Highs_call, "libhighs.so"), Cint, (Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ref{Cint}),
   n_col, n_row, n_nz, colcost, collower, colupper, rowlower, rowupper, matstart, matindex, matvalue, solution.colvalue, solution.coldual, solution.rowvalue, solution.rowdual, basis.colbasisstatus, basis.rowbasisstatus, modelstatus)

   return status, solution, basis, modelstatus
end

function Highs_create()
   return ccall((:Highs_create, "libhighs.so"), Ptr{Cvoid}, () )
end

function Highs_destroy(highs)
   ccall((:Highs_destroy, "libhighs.so"), Cvoid, (Ptr{Cvoid},), highs)
end

function Highs_run(highs)
   return ccall((:Highs_run, "libhighs.so"), Cint, (Ptr{Cvoid},), highs)
end

function Highs_readModel(highs, filename)
   name = convert(Cstring, pointer(filename))
   return ccall((:Highs_readModel, "libhighs.so"), Cint, (Ptr{Cvoid}, Cstring), highs, name)
end

function Highs_writeModel(highs, filename)
   name = convert(Cstring, pointer(filename))
   return ccall((:Highs_writeModel, "libhighs.so"), Cint, (Ptr{Cvoid}, Cstring), highs, name)
end

function Highs_passLp(highs, model)
   n_col = convert(Cint, size(model.colcost, 1))
   n_row = convert(Cint, size(model.rowlower, 1))
   n_nz = convert(Cint, size(model.aindex, 1))

   colcost = convert(Array{Cdouble}, model.colcost)
   collower = convert(Array{Cdouble}, model.collower)
   colupper = convert(Array{Cdouble}, model.colupper)

   rowlower = convert(Array{Cdouble}, model.rowlower)
   rowupper = convert(Array{Cdouble}, model.rowupper)
   matstart = convert(Array{Cint}, model.astart)
   matindex = convert(Array{Cint}, model.aindex)
   matvalue = convert(Array{Cdouble}, model.avalue)

   return ccall((:Highs_passLp, "libhighs.so"), Cint, (Ptr{Cvoid},Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint},Ptr{Cint},Ptr{Cdouble}),
   highs, n_col, n_row, n_nz, colcost, collower, colupper, rowlower, rowupper, matstart, matindex, matvalue)
end

function Highs_setOptionValue(highs, option, value)
   opt = convert(Cstring, pointer(option))
   val = convert(Cstring, pointer(val))
   return ccall((:Highs_setOptionValue, "libhighs.so"), Cint, (Ptr{Cvoid}, Cstring, Cstring), highs, opt, val)
end

function Highs_getNumCols(highs)
   return ccall((:Highs_getNumCols, "libhighs.so"), Cint, ())
end

function Highs_getNumRows(highs)
   return ccall((:Highs_getNumRows, "libhighs.so"), Cint, ())
end

function Highs_getNumNz(highs)
   return ccall((:Highs_getNumNz, "libhighs.so"), Cint, ())
end

function Highs_getSolution(highs)
   numcol = Highs_getNumCols(highs)
   numrow = Highs_getNumRows(highs)

   solution = HighsSolution(Array{Cdouble, 1}(undef, numcol), Array{Cdouble, 1}(undef, numcol), Array{Cdouble, 1}(undef, numrow),  Array{Cdouble, 1}(undef, numrow))

   ccall((:Highs_getSolution, "libhighs.so"), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, solution.colvalue, solution.coldual, solution.rowvalue, solution.rowdual)

   return solution
end

function Highs_getBasis(highs)
   numcol = Highs_getNumCols(highs)
   numrow = Highs_getNumRows(highs)

   basis = HighsBasis(Array{Cint, 1}(undef, numcol), Array{Cint, 1}(undef, numrow))   

   ccall((:Highs_getBasis, "libhighs.so"), Cvoid, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}), highs, basis.colbasisstatus, basis.rowbasisstatus)

   return basis
end

function Highs_getObjectiveValue(highs)
   return ccall((:Highs_getObjectiveValue, "libhighs.so"), Cdouble, (Ptr{Cvoid},), highs)
end

function Highs_getIterationCount(highs)
   return ccall((:Highs_getIterationCount, "libhighs.so"), Cint, (Ptr{Cvoid},), highs)
end

function Highs_addRow(highs, lower, upper, indices, values)
   n_new_nz = convert(Cint, size(indices, 1))

   lo = convert(Cdouble, lower)
   hi = convert(Cdouble, upper)

   idx = convert(Array{Cint}, indices)
   val = convert(Array{Cdouble}, values)

   return ccall((:Highs_addRow, "libhighs.so"), Cint, (Ptr{Cvoid}, Cdouble, Cdouble, Cint, Ptr{Cint}, Ptr{Cdouble}), highs, lo, hi, n_new_nz, idx, val)
end

function Highs_addRows(highs, lower, upper, starts, indices, values)
   n_new_rows = convert(Cint, size(lower, 1))
   n_new_nz = convert(Cint, size(indices, 1))

   lo = convert(Array{Cdouble}, lower)
   hi = convert(Array{Cdouble}, upper)

   st = convert(Array{Cint}, starts)
   idx = convert(Array{Cint}, indices)
   val = convert(Array{Cdouble}, values)

   return ccall((:Highs_addRows, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), 
   highs, n_new_rows, lo, hi, n_new_nz, st, idx, val)
end

function Highs_addCol(highs, cost, lower, upper, indices, values)
   n_new_nz = convert(Cint, size(indices, 1))

   cc = convert(Cdouble, cost)
   lo = convert(Cdouble, lower)
   hi = convert(Cdouble, upper)

   idx = convert(Array{Cint}, indices)
   val = convert(Array{Cdouble}, values)

   return ccall((:Highs_addCol, "libhighs.so"), Cint, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cint, Ptr{Cint}, Ptr{Cdouble}),
   highs, cc, lo, hi, n_new_nz, idx, val)
end

function Highs_addCols(highs, costs, lower, upper, starts, indices, values)
   n_new_cols = convert(Cint, size(lower, 1))
   n_new_nz = convert(Cint, size(indices, 1))

   co = convert(Array{Cdouble}, costs)
   lo = convert(Array{Cdouble}, lower)
   hi = convert(Array{Cdouble}, upper)

   st = convert(Array{Cint}, starts)
   idx = convert(Array{Cint}, indices)
   val = convert(Array{Cdouble}, values)

   return ccall((:Highs_addCols, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), 
   highs, n_new_cols, co, lo, hi, n_new_nz, st, idx, val)
end

function Highs_changeObjectiveSense(highs, sense)
   sns = convert(Cint, sense)

   return ccall((:Highs_changeObjectiveSense, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint), highs, sns)
end

function Highs_changeColCost(highs, colidx, cost)
   col = convert(Cint, colidx)
   cst = convert(Cdouble, cost)

   return ccall((:Highs_changeColCost, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Cdouble), highs, colidx, cost)
end

function Highs_changeColsCostBySet(highs, set, cost)
   num_set_entries = convert(Cint, size(set, 1)) 

   st = convert(Array{Cint}, set)
   cst = convert(Array{Cdouble}, cost)

   return ccall((:Highs_changeColsCostBySet, "libhighs.so"), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cdouble}), highs, st, cst)
end

function Highs_changeColsCostByMask(highs, mask, cost)
   msk = convert(Array{Cint}, mask)
   cst = convert(Array{Cdouble}, cost)

   return ccall((:Highs_changeColsCostByMask, "libhighs.so"), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cdouble}), highs, msk, cst)
end

function Highs_changeColBounds(highs, col, lower, upper)
   colidx = convert(Cint, col)

   lo = convert(Cdouble, lower)
   hi = convert(Cdouble, upper)

   return ccall((:Highs_changeColBounds, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Cdouble, Cdouble), highs, colidx, lo, hi)
end

function Highs_changeColsBoundsByRange(highs, from, to, lower, upper)
   f = convert(Cint, from)
   t = convert(Int33, to)

   lo = convert(Array{Cdouble}, lower)
   hi = convert(Array{Cdouble}, upper)

   return ccall((:Highs_changeColsBoundsByRange, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
   highs, f, t, lo, hi)
end

function Highs_changeColsBoundsBySet(highs, set, lower, upper)
   nset = convert(Cint, size(set, 1))

   st = convert(Array{Cint}, set)

   lo = convert(Array{Cdouble}, lower)
   hi = convert(Array{Cdouble}, upper)

   return ccall((:Highs_changeColsBoundsBySet, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), 
   highs, nset, st, lo, hi)
end

function Highs_changeColsBoundsByMask(highs, mask, lower, upper) 
   msk = convert(Array{Cint}, mask)

   lo = convert(Array{Cdouble}, lower)
   hi = convert(Array{Cdouble}, upper)

   return ccall((:Highs_changeColsBoundsByMask, "libhighs.so"), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}),
   highs, msk, lo, hi)
end

function Highs_changeRowBounds(highs, row, lower, upper) 
   idx = convert(Cint, row)
   lo = convert(Cdouble, lower)
   hi = convert(Cdouble, upper)

   return ccall((:Highs_changeRowBounds, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Cdouble, Cdouble),
   highs, idx, lo, hi)
end

function Highs_changeRowsBoundsBySet(highs, set, lower, upper)
   nst = convert(Cint, size(set, 1))
   st = convert(Array{Cint}, set)

   lo = convert(Array{Cdouble}, lower)
   hi = convert(Array(Cdouble), upper)

   return ccall((:Highs_changeRowsBoundsBySet, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}),
   highs, nst, st, lo, hi)
end

function Highs_changeRowsBoundsByMask(highs, mask, lower, upper)
   msk = convert(Array{Cint}, mask)

   lo = convert(Array{Cdouble}, lower)
   hi = convert(Array{Cdouble}, upper)

   return ccall((:Highs_changeRowsBoundsByMask, "libhighs.so"), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}),
   highs, msk, lo, hi)
end

#=

int Highs_getColsByRange(
    void *highs,          //!< HiGHS object reference
    const int from_col,   //!< The index of the first column to
                          //!< get from the model
    const int to_col,     //!< One more than the last column to get
                          //!< from the model
    int* num_col,          //!< Number of columns got from the model
    double *costs,        //!< Array of size num_col with costs
    double *lower,        //!< Array of size num_col with lower bounds
    double *upper,        //!< Array of size num_col with upper bounds
    int* num_nz,           //!< Number of nonzeros got from the model
    int *matrix_start,    //!< Array of size num_col with start
                          //!< indices of the columns
    int *matrix_index,    //!< Array of size num_nz with row
                          //!< indices for the columns
    double *matrix_value  //!< Array of size num_nz with row
                          //!< values for the columns
);

int Highs_getColsBySet(
    void *highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int *set,             //!< Array of size num_set_entries with indices
                                //!< of columns to get
    int* num_col,                //!< Number of columns got from the model
    double *costs,              //!< Array of size num_col with costs
    double *lower,              //!< Array of size num_col with lower bounds
    double *upper,              //!< Array of size num_col with upper bounds
    int* num_nz,                 //!< Number of nonzeros got from the model
    int *matrix_start,          //!< Array of size num_col with start indices
                                //!< of the columns
    int *matrix_index,          //!< Array of size num_nz with row indices
                                //!< for the columns
    double *matrix_value        //!< Array of size num_nz with row values
                                //!< for the columns
);

int Highs_getColsByMask(
    void *highs,          //!< HiGHS object reference
    const int *mask,      //!< Full length array with 1 => get; 0 => not
    int* num_col,          //!< Number of columns got from the model
    double *costs,        //!< Array of size num_col with costs
    double *lower,        //!< Array of size num_col with lower bounds
    double *upper,        //!< Array of size num_col with upper bounds
    int* num_nz,           //!< Number of nonzeros got from the model
    int *matrix_start,    //!<  Array of size num_col with start
                          //!<  indices of the columns
    int *matrix_index,    //!<  Array of size num_nz with row indices
                          //!<  for the columns
    double *matrix_value  //!<  Array of size num_nz with row values
                          //!<  for the columns
);

int Highs_getRowsByRange(
    void *highs,          //!< HiGHS object reference
    const int from_row,   //!< The index of the first row to get from the model
    const int to_row,     //!< One more than the last row get from the model
    int* num_row,          //!< Number of rows got from the model
    double *lower,        //!< Array of size num_row with lower bounds
    double *upper,        //!< Array of size num_row with upper bounds
    int* num_nz,           //!< Number of nonzeros got from the model
    int *matrix_start,    //!< Array of size num_row with start indices of the
                          //!< rows
    int *matrix_index,    //!< Array of size num_nz with column indices for the
                          //!< rows
    double *matrix_value  //!< Array of size num_nz with column values for the
                          //!< rows
);

int Highs_getRowsBySet(
    void *highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int *set,             //!< Array of size num_set_entries with indices
                                //!< of rows to get
    int* num_row,                //!< Number of rows got from the model
    double *lower,              //!< Array of size num_row with lower bounds
    double *upper,              //!< Array of size num_row with upper bounds
    int* num_nz,                 //!< Number of nonzeros got from the model
    int *matrix_start,          //!< Array of size num_row with start indices
                                //!< of the rows
    int *matrix_index,          //!< Array of size num_nz with column indices
                                //!< for the rows
    double *matrix_value        //!< Array of size num_nz with column
                                //!< values for the rows
);

int Highs_getRowsByMask(
    void *highs,          //!< HiGHS object reference
    const int *mask,      //!< Full length array with 1 => get; 0 => not
    int* num_row,          //!< Number of rows got from the model
    double *lower,        //!< Array of size num_row with lower bounds
    double *upper,        //!< Array of size num_row with upper bounds
    int* num_nz,           //!< Number of nonzeros got from the model
    int *matrix_start,    //!< Array of size num_row with start indices
                          //!< of the rows
    int *matrix_index,    //!< Array of size num_nz with column indices
                          //!< for the rows
    double *matrix_value  //!< Array of size num_nz with column
                          //!< values for the rows
);

=#

function Highs_deleteColsByRange(highs, from, to)
   f = convert(Cint, from)
   t = convert(Cint, to)

   return ccall((:Highs_deleteColsByRange, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Cint), highs, f, t)
end

function Highs_deleteColsBySet(highs, set)
   nset = convert(Cint, size(set, 1))

   st = convert(Array{Cint}, set)

   return ccall((:Highs_deleteColsBySet, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}), highs, nset, st)
end

function Highs_deleteColsByMask(highs, mask)
   msk = convert(Array{Cint}, mask)

   return ccall((:Highs_deleteColsByMask, "libighs.so"), Cint, (Ptr{Cvoid}, Ptr{Cint}), highs, msk)
end

function Highs_deleteRowsByRange(highs, from, to)
   f = convert(Cint, from)
   t = convert(Cint, to)   

   return ccall((:Highs_deleteRowsByRange, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Cint), highs, f, t)
end

function Highs_deleteRowsBySet(highs, set)
   nset = convert(Cint, size(set, 1))

   st = convert(Array{Cint}, set)

   return ccall((:Highs_deleteRowsBySet, "libhighs.so"), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}), highs, nset, st)
end

function Highs_deleteRowsByMask(highs, mask)
   msk = convert(Array{Cint}, mask)

   return ccall((:Highs_deleteRowsByMask, "libhighs.so"), Cint, (Ptr{Cvoid}, Ptr{Cint}), highs, msk)
end
