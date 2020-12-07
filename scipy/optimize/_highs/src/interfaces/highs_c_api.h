/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_C_API
#define HIGHS_C_API

#ifdef __cplusplus
extern "C" {
#endif

/*
 * @brief runs a model using HiGHS
 */
int Highs_call(
    int numcol,        //!< number of columns
    int numrow,        //!< number of rows
    int numnz,         //!< number of entries in the constraint matrix
    double* colcost,   //!< array of length [numcol] with column costs
    double* collower,  //!< array of length [numcol] with lower column bounds
    double* colupper,  //!< array of length [numcol] with upper column bounds
    double* rowlower,  //!< array of length [numrow] with lower row bounds
    double* rowupper,  //!< array of length [numrow] with upper row bounds
    int* astart,       //!< array of length [numcol+1] with column start indices
    int*
        aindex,  //!< array of length [numnz] with row indices of matrix entries
    double* avalue,    //!< array of length [numnz] with value of matrix entries
    double* colvalue,  //!< array of length [numcol], filled with column values
    double* coldual,   //!< array of length [numcol], filled with column duals
    double* rowvalue,  //!< array of length [numrow], filled with row values
    double* rowdual,   //!< array of length [numrow], filled with row duals
    int* colbasisstatus,  //!< array of length [numcol], filled with column
                          //!< basis stati
    int* rowbasisstatus,  //!< array of length [numrow], filled with row basis
                          //!< stati
    int* modelstatus      //!< status of the model will be saved here
);

/*
 * @brief creates a HiGHS object and returns the reference
 */
void* Highs_create(void);

/*
 * @brief destroys a HiGHS object
 */
void Highs_destroy(void* highs  //!< HiGHS object reference
);

/*
 * @brief
 */
int Highs_readModel(void* highs,          //!< HiGHS object reference
                    const char* filename  //!< filename
);

/*
 * @brief
 */
int Highs_writeModel(void* highs,          //!< HiGHS object reference
                     const char* filename  //!< filename
);

/*
 * @brief
 */
int Highs_clearModel(void* highs  //!< HiGHS object reference
);

/*
 * @brief
 */
int Highs_runQuiet(void* highs  //!< HiGHS object reference
);

/*
 * @brief Runs HiGHS
 */
int Highs_run(void* highs  //!< HiGHS object reference
);

/*
 * @brief Reports the solution and basis status
 */
int Highs_writeSolution(void* highs,          //!< HiGHS object reference
                        const char* filename  //!< filename
);

/*
 * @brief pass an LP to HiGHS
 */
int Highs_passLp(
    void* highs,       //!< HiGHS object reference
    int numcol,        //!< number of columns
    int numrow,        //!< number of rows
    int numnz,         //!< number of entries in the constraint matrix
    double* colcost,   //!< array of length [numcol] with column costs
    double* collower,  //!< array of length [numcol] with lower column bounds
    double* colupper,  //!< array of length [numcol] with upper column bounds
    double* rowlower,  //!< array of length [numrow] with lower row bounds
    double* rowupper,  //!< array of length [numrow] with upper row bounds
    int* astart,       //!< array of length [numcol+1] with column start indices
    int*
        aindex,  //!< array of length [numnz] with row indices of matrix entries
    double* avalue  //!< array of length [numnz] with value of matrix entries
);

/*
 * @brief
 */
int Highs_setHighsBoolOptionValue(void* highs,  //!< HiGHS object reference
                                  const char* option,  //!< name of the option
                                  const int value      //!< new value of option
);

/*
 * @brief
 */
int Highs_setHighsIntOptionValue(void* highs,  //!< HiGHS object reference
                                 const char* option,  //!< name of the option
                                 const int value      //!< new value of option
);

/*
 * @brief
 */
int Highs_setHighsDoubleOptionValue(void* highs,  //!< HiGHS object reference
                                    const char* option,  //!< name of the option
                                    const double value  //!< new value of option
);

/*
 * @brief
 */
int Highs_setHighsStringOptionValue(void* highs,  //!< HiGHS object reference
                                    const char* option,  //!< name of the option
                                    const char* value  //!< new value of option
);

/*
 * @brief
 */
int Highs_setHighsOptionValue(void* highs,         //!< HiGHS object reference
                              const char* option,  //!< name of the option
                              const char* value    //!< new value of option
);

/*
 * @brief
 */
int Highs_getHighsBoolOptionValue(void* highs,  //!< HiGHS object reference
                                  const char* option,  //!< name of the option
                                  int* value           //!< value of option
);

/*
 * @brief
 */
int Highs_getHighsIntOptionValue(void* highs,  //!< HiGHS object reference
                                 const char* option,  //!< name of the option
                                 int* value           //!< value of option
);

/*
 * @brief
 */
int Highs_getHighsDoubleOptionValue(void* highs,  //!< HiGHS object reference
                                    const char* option,  //!< name of the option
                                    double* value        //!< value of option
);

/*
 * @brief
 */
int Highs_getHighsStringOptionValue(
    void* highs,         //!< HiGHS object reference
    const char* option,  //!< name of the option
    char* value  //!< pointer to allocated memory to store value of option
);

/*
 * @brief
 */
int Highs_resetHighsOptions(void* highs  //!< HiGHS object reference
);

/*
 * @brief
 */
int Highs_getHighsIntInfoValue(void* highs,       //!< HiGHS object reference
                               const char* info,  //!< The info name
                               int* value         //!< The info value
);

/*
 * @brief
 */
int Highs_getHighsDoubleInfoValue(void* highs,       //!< HiGHS object reference
                                  const char* info,  //!< The info name
                                  double* value      //!< The info value
);
/*
 * @brief
 */
void Highs_getSolution(
    void* highs,       //!< HiGHS object reference
    double* colvalue,  //!< array of length [numcol], filled with column values
    double* coldual,   //!< array of length [numcol], filled with column duals
    double* rowvalue,  //!< array of length [numrow], filled with row values
    double* rowdual    //!< array of length [numrow], filled with row duals
);

/*
 * @brief
 */
void Highs_getBasis(
    void* highs,     //!< HiGHS object reference
    int* colstatus,  //!< array of length [numcol], filled with column basis
                     //!< stati
    int* rowstatus   //!< array of length [numrow], filled with row basis stati
);

/**
 * @brief Returns the status of the (scaled) model
 */
int Highs_getModelStatus(
    void* highs,            //!< HiGHS object reference
    const int scaled_model  //!< 0 (nonzero) for status of (scaled) model
);

/**
 * @brief Returns the objective function value (if known)
 */
double Highs_getObjectiveValue(void* highs  //!< HiGHS object reference,
);

/**
 * @brief Returns the simplex iteration count (if known)
 */
int Highs_getIterationCount(void* highs  //!< HiGHS object reference,
);

/**
 * @brief Returns the simplex iteration count (if known)
 */
int Highs_getSimplexIterationCount(void* highs  //!< HiGHS object reference,
);

/**
 * @brief Gets the basic variables in the order corresponding to
 * calls to getBasisInverseRow, getBasisInverseCol, getBasisSolve,
 * getBasisTransposeSolve, getReducedRow and getReducedColumn. As
 * required by SCIP, non-negative entries are indices of columns,
 * and negative entries are -(row_index+1).
 */
int Highs_getBasicVariables(void* highs,          //!< HiGHS object reference,
                            int* basic_variables  //!< Basic variables
);

/**
 * @brief Gets a row of \f$B^{-1}\f$ for basis matrix \f$B\f$
 */
int Highs_getBasisInverseRow(void* highs,         //!< HiGHS object reference
                             const int row,       //!< Index of row required
                             double* row_vector,  //!< Row required
                             int* row_num_nz,     //!< Number of nonzeros
                             int* row_indices     //!< Indices of nonzeros
);

/**
 * @brief Gets a column of \f$B^{-1}\f$ for basis matrix \f$B\f$
 */
int Highs_getBasisInverseCol(void* highs,         //!< HiGHS object reference
                             const int col,       //!< Index of column required
                             double* col_vector,  //!< Column required
                             int* col_num_nz,     //!< Number of nonzeros
                             int* col_indices     //!< Indices of nonzeros
);

/**
 * @brief Forms \f$\mathbf{x}=B^{-1}\mathbf{b}\f$ for a given vector
 * \f$\mathbf{b}\f$
 */
int Highs_getBasisSolve(void* highs,              //!< HiGHS object reference
                        const double* rhs,        //!< RHS \f$\mathbf{b}\f$
                        double* solution_vector,  //!< Solution \f$\mathbf{x}\f$
                        int* solution_num_nz,     //!< Number of nonzeros
                        int* solution_indices     //!< Indices of nonzeros
);

/**
 * @brief Forms \f$\mathbf{x}=B^{-T}\mathbf{b}\f$ for a given vector
 * \f$\mathbf{b}\f$
 */
int Highs_getBasisTransposeSolve(
    void* highs,              //!< HiGHS object reference
    const double* rhs,        //!< RHS \f$\mathbf{b}\f$
    double* solution_vector,  //!< Solution  \f$\mathbf{x}\f$
    int* solution_nz,         //!< Number of nonzeros
    int* solution_indices     //!< Indices of nonzeros
);

/**
 * @brief Forms a row of \f$B^{-1}A\f$
 */
int Highs_getReducedRow(void* highs,         //!< HiGHS object reference
                        const int row,       //!< Index of row required
                        double* row_vector,  //!< Row required
                        int* row_num_nz,     //!< Number of nonzeros
                        int* row_indices     //!< Indices of nonzeros
);

/**
 * @brief Forms a column of \f$B^{-1}A\f$
 */
int Highs_getReducedColumn(void* highs,         //!< HiGHS object reference
                           const int col,       //!< Index of column required
                           double* col_vector,  //!< Column required
                           int* col_num_nz,     //!< Number of nonzeros
                           int* col_indices     //!< Indices of nonzeros
);

/**
 * @brief Passes a basis to HiGHS
 */
int Highs_setBasis(void* highs,           //!< HiGHS object reference
                   const int* colstatus,  //!< Column status
                   const int* rowstatus   //!< Row status
);

/**
 * @brief Sets up a logical basis in HiGHS
 */
int Highs_setLogicalBasis(void* highs  //!< HiGHS object reference
);

/**
 * @brief Returns the cumulative wall-clock time spent in Highs_run();
 */
double Highs_getHighsRunTime(void* highs  //!< HiGHS object reference
);

/**
 * @brief Adds a row to the model
 */
int Highs_addRow(
    void* highs,           //!< HiGHS object reference
    const double lower,    //!< Lower bound of the row
    const double upper,    //!< Upper bound of the row
    const int num_new_nz,  //!< Number of nonzeros in the row
    const int* indices,    //!< Array of size num_new_nz with column indices
    const double* values   //!< Array of size num_new_nz with column values
);

/**
 * @brief Adds multiple rows to the model
 */
int Highs_addRows(
    void* highs,            //!< HiGHS object reference
    const int num_new_row,  //!< Number of new rows
    const double* lower,    //!< Array of size num_new_row with lower bounds
    const double* upper,    //!< Array of size num_new_row with upper bounds
    const int num_new_nz,   //!< Number of new nonzeros
    const int*
        starts,  //!< Array of size num_new_row with start indices of the rows
    const int*
        indices,  //!< Array of size num_new_nz with column indices for all rows
    const double*
        values  //!< Array of size num_new_nz with column values for all rows
);

/**
 * @brief Adds a column to the model
 */
int Highs_addCol(
    void* highs,           //!< HiGHS object reference
    const double cost,     //!< Cost of the column
    const double lower,    //!< Lower bound of the column
    const double upper,    //!< Upper bound of the column
    const int num_new_nz,  //!< Number of nonzeros in the column
    const int* indices,    //!< Array of size num_new_nz with row indices
    const double* values   //!< Array of size num_new_nz with row values
);

/**
 * @brief Adds multiple columns to the model
 */
int Highs_addCols(
    void* highs,            //!< HiGHS object reference
    const int num_new_col,  //!< Number of new columns
    const double* costs,    //!< Array of size num_new_col with costs
    const double* lower,    //!< Array of size num_new_col with lower bounds
    const double* upper,    //!< Array of size num_new_col with upper bounds
    const int num_new_nz,   //!< Number of new nonzeros
    const int* starts,      //!< Array of size num_new_row with start indices of
                            //!< the columns
    const int* indices,     //!< Array of size num_new_nz with row indices for
                            //!< all columns
    const double*
        values  //!< Array of size num_new_nz with row values for all columns
);

/**
 * @brief Change the objective sense of the model
 */
int Highs_changeObjectiveSense(void* highs,     //!< HiGHS object reference
                               const int sense  //!< New objective sense
);

/**
 * @brief Change the cost of a column
 */
int Highs_changeColCost(
    void* highs,       //!< HiGHS object reference
    const int col,     //!< The index of the column whose cost is to change
    const double cost  //!< The new cost
);

/**
 * @brief Change the cost of multiple columns given by a set of indices
 */
int Highs_changeColsCostBySet(
    void* highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int* set,     //!< Array of size num_set_entries with indices of
                        //!< columns whose costs change
    const double* cost  //!< Array of size num_set_entries with new costs
);

/**
 * @brief Change the cost of multiple columns given by a mask
 */
int Highs_changeColsCostByMask(
    void* highs,        //!< HiGHS object reference
    const int* mask,    //!< Full length array with 1 => change; 0 => not
    const double* cost  //!< Full length array of new costs
);

/**
 * @brief Change the bounds of a column
 */
int Highs_changeColBounds(
    void* highs,         //!< HiGHS object reference
    const int col,       //!< The index of the column whose bounds are to change
    const double lower,  //!< The new lower bound
    const double upper   //!< The new upper bound
);

/**
 * @brief Change the bounds of multiple columns given by an interval
 */
int Highs_changeColsBoundsByRange(
    void* highs,         //!< HiGHS object reference
    const int from_col,  //!< The index of the first column whose bounds change
    const int to_col,    //!< One more than the index of the last column whose
                         //!< bounds change
    const double*
        lower,  //!< Array of size to_col-from_col with new lower bounds
    const double*
        upper  //!< Array of size to_col-from_col with new upper bounds
);

/**
 * @brief Change the bounds of multiple columns given by a set of indices
 */
int Highs_changeColsBoundsBySet(
    void* highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int* set,  //!< Array of size num_set_entries with indices of
                     //!< columns whose bounds change
    const double*
        lower,  //!< Array of size num_set_entries with new lower bounds
    const double*
        upper  //!< Array of size num_set_entries with new upper bounds
);

/**
 * @brief Change the cost of multiple columns given by a mask
 */
int Highs_changeColsBoundsByMask(
    void* highs,          //!< HiGHS object reference
    const int* mask,      //!< Full length array with 1 => change; 0 => not
    const double* lower,  //!< Full length array of new lower bounds
    const double* upper   //!< Full length array of new upper bounds
);

/**
 * @brief Change the bounds of a row
 */
int Highs_changeRowBounds(
    void* highs,         //!< HiGHS object reference
    const int row,       //!< The index of the row whose bounds are to change
    const double lower,  //!< The new lower bound
    const double upper   //!< The new upper bound
);

/**
 * @brief Change the bounds of multiple rows given by a set of indices
 */
int Highs_changeRowsBoundsBySet(
    void* highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int* set,  //!< Array of size num_set_entries with indices of rows
                     //!< whose bounds change
    const double*
        lower,  //!< Array of size num_set_entries with new lower bounds
    const double*
        upper  //!< Array of size num_set_entries with new upper bounds
);

/**
 * @brief Change the cost of multiple rows given by a mask
 */
int Highs_changeRowsBoundsByMask(
    void* highs,          //!< HiGHS object reference
    const int* mask,      //!< Full length array with 1 => change; 0 => not
    const double* lower,  //!< Full length array of new lower bounds
    const double* upper   //!< Full length array of new upper bounds
);

/**
 * @brief Get the objective sense
 */
int Highs_getObjectiveSense(void* highs, int* sense);
/**
 * @brief Get multiple columns from the model given by an interval
 */
int Highs_getColsByRange(
    void* highs,          //!< HiGHS object reference
    const int from_col,   //!< The index of the first column to
                          //!< get from the model
    const int to_col,     //!< One more than the last column to get
                          //!< from the model
    int* num_col,         //!< Number of columns got from the model
    double* costs,        //!< Array of size num_col with costs
    double* lower,        //!< Array of size num_col with lower bounds
    double* upper,        //!< Array of size num_col with upper bounds
    int* num_nz,          //!< Number of nonzeros got from the model
    int* matrix_start,    //!< Array of size num_col with start
                          //!< indices of the columns
    int* matrix_index,    //!< Array of size num_nz with row
                          //!< indices for the columns
    double* matrix_value  //!< Array of size num_nz with row
                          //!< values for the columns
);

/**
 * @brief Get multiple columns from the model given by a set
 */
int Highs_getColsBySet(
    void* highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int* set,             //!< Array of size num_set_entries with indices
                                //!< of columns to get
    int* num_col,               //!< Number of columns got from the model
    double* costs,              //!< Array of size num_col with costs
    double* lower,              //!< Array of size num_col with lower bounds
    double* upper,              //!< Array of size num_col with upper bounds
    int* num_nz,                //!< Number of nonzeros got from the model
    int* matrix_start,          //!< Array of size num_col with start indices
                                //!< of the columns
    int* matrix_index,          //!< Array of size num_nz with row indices
                                //!< for the columns
    double* matrix_value        //!< Array of size num_nz with row values
                                //!< for the columns
);

/**
 * @brief Get multiple columns from the model given by a mask
 */
int Highs_getColsByMask(
    void* highs,          //!< HiGHS object reference
    const int* mask,      //!< Full length array with 1 => get; 0 => not
    int* num_col,         //!< Number of columns got from the model
    double* costs,        //!< Array of size num_col with costs
    double* lower,        //!< Array of size num_col with lower bounds
    double* upper,        //!< Array of size num_col with upper bounds
    int* num_nz,          //!< Number of nonzeros got from the model
    int* matrix_start,    //!<  Array of size num_col with start
                          //!<  indices of the columns
    int* matrix_index,    //!<  Array of size num_nz with row indices
                          //!<  for the columns
    double* matrix_value  //!<  Array of size num_nz with row values
                          //!<  for the columns
);

/**
 * @brief Get multiple rows from the model given by an interval
 */
int Highs_getRowsByRange(
    void* highs,          //!< HiGHS object reference
    const int from_row,   //!< The index of the first row to get from the model
    const int to_row,     //!< One more than the last row get from the model
    int* num_row,         //!< Number of rows got from the model
    double* lower,        //!< Array of size num_row with lower bounds
    double* upper,        //!< Array of size num_row with upper bounds
    int* num_nz,          //!< Number of nonzeros got from the model
    int* matrix_start,    //!< Array of size num_row with start indices of the
                          //!< rows
    int* matrix_index,    //!< Array of size num_nz with column indices for the
                          //!< rows
    double* matrix_value  //!< Array of size num_nz with column values for the
                          //!< rows
);

/**
 * @brief Get multiple rows from the model given by a set
 */
int Highs_getRowsBySet(
    void* highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int* set,             //!< Array of size num_set_entries with indices
                                //!< of rows to get
    int* num_row,               //!< Number of rows got from the model
    double* lower,              //!< Array of size num_row with lower bounds
    double* upper,              //!< Array of size num_row with upper bounds
    int* num_nz,                //!< Number of nonzeros got from the model
    int* matrix_start,          //!< Array of size num_row with start indices
                                //!< of the rows
    int* matrix_index,          //!< Array of size num_nz with column indices
                                //!< for the rows
    double* matrix_value        //!< Array of size num_nz with column
                                //!< values for the rows
);

/**
 * @brief Get multiple rows from the model given by a mask
 */
int Highs_getRowsByMask(
    void* highs,          //!< HiGHS object reference
    const int* mask,      //!< Full length array with 1 => get; 0 => not
    int* num_row,         //!< Number of rows got from the model
    double* lower,        //!< Array of size num_row with lower bounds
    double* upper,        //!< Array of size num_row with upper bounds
    int* num_nz,          //!< Number of nonzeros got from the model
    int* matrix_start,    //!< Array of size num_row with start indices
                          //!< of the rows
    int* matrix_index,    //!< Array of size num_nz with column indices
                          //!< for the rows
    double* matrix_value  //!< Array of size num_nz with column
                          //!< values for the rows
);

/**
 * @brief Delete multiple columns from the model given by an interval
 */
int Highs_deleteColsByRange(
    void* highs,         //!< HiGHS object reference
    const int from_col,  //!< The index of the first column
                         //!< to delete from the model
    const int to_col     //!< One more than the last column to
                         //!< delete from the model
);

/**
 * @brief Delete multiple columns from the model given by a set
 */
int Highs_deleteColsBySet(
    void* highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int* set  //!< Array of size num_set_entries with indices of columns
                    //!< to delete
);

/**
 * @brief Delete multiple columns from the model given by a mask
 */
int Highs_deleteColsByMask(
    void* highs,  //!< HiGHS object reference
    int* mask     //!< Full length array with 1 => delete; 0 => not
);

/**
 * @brief Delete multiple rows from the model given by an interval
 */
int Highs_deleteRowsByRange(
    void* highs,  //!< HiGHS object reference
    const int
        from_row,     //!< The index of the first row to delete from the model
    const int to_row  //!< One more than the last row delete from the model
);

/**
 * @brief Delete multiple rows from the model given by a set
 */
int Highs_deleteRowsBySet(
    void* highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int* set  //!< Array of size num_set_entries with indices of columns
                    //!< to delete
);

/**
 * @brief Delete multiple rows from the model given by a mask
 */
int Highs_deleteRowsByMask(
    void* highs,  //!< HiGHS object reference
    int* mask     //!< Full length array with 1 => delete; 0 => not
);

/**
 * @brief Returns the value of infinity used by HiGHS
 */
double Highs_getHighsInfinity(void* highs  //!< HiGHS object reference
);

/**
 * @brief Returns the number of columns of the current model
 */
int Highs_getNumCols(void* highs  //!< HiGHS object reference
);

/**
 * @brief Returns the number of rows of the current model
 */
int Highs_getNumRows(void* highs  //!< HiGHS object reference
);

/**
 * @brief Returns the number of nonzeroes of the current model
 */
int Highs_getNumNz(void* highs  //!< HiGHS object reference
);

/**
 * @brief Returns a pointer to a character representation of a HiGHS model
 * status
 */
const char* Highs_highsModelStatusToChar(
    void* highs,                //!< HiGHS object reference
    int int_highs_model_status  //!< Status to interpret
);

/**
 * @brief Returns a pointer to a character representation of a primal/dual
 * status
 */
const char* Highs_primalDualStatusToChar(
    void* highs,                //!< HiGHS object reference
    int int_primal_dual_status  //!< Status to interpret
);

// /**
//  * @brief Returns the current model
//  */
// void Highs_getLp(
//     void *highs,       //!< HiGHS object reference
//     int* numcol,        //!< number of columns
//     int* numrow,        //!< number of rows
//     int* numnz,         //!< number of entries in the constraint matrix
//     double *colcost,   //!< array of length [numcol] with column costs
//     double *collower,  //!< array of length [numcol] with lower column bounds
//     double *colupper,  //!< array of length [numcol] with upper column bounds
//     double *rowlower,  //!< array of length [numrow] with lower row bounds
//     double *rowupper,  //!< array of length [numrow] with upper row bounds
//     int *astart,       //!< array of length [numcol+1] with column start
//     indices int *
//         aindex,  //!< array of length [numnz] with row indices of matrix
//         entries
//     double *avalue  //!< array of length [numnz] with value of matrix entries
// );

// Highs(HighsOptions &options) { options_ = options; }

// /**
//  * @brief Clears the vector of HighsModelObjects (hmos), creates a
//  * HighsModelObject for this LP and makes it the first of the vector
//  * of HighsModelObjects
//  */
// HighsStatus passModel(
//     const HighsLp &lp  //!< The HighsLp instance for this LP
// );

// HighsStatus setSolution(const HighsSolution &solution);

// /**
//  * @brief Uses the HighsBasis passed to set the basis for the
//  * LP of the (first?) HighsModelObject
//  */
// HighsStatus setBasis(const HighsBasis &basis);

#ifdef __cplusplus
}
#endif

#endif
