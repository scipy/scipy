/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HMPSIO.h
 * @brief
 */
#ifndef IO_HMPSIO_H_
#define IO_HMPSIO_H_

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "io/Filereader.h"
#include "util/HighsType.h"

using std::string;
using std::vector;

const HighsInt MPS_ROW_TY_N = 0;
const HighsInt MPS_ROW_TY_E = 1;
const HighsInt MPS_ROW_TY_L = 2;
const HighsInt MPS_ROW_TY_G = 3;
const HighsInt field_1_start = 1;
const HighsInt field_1_width = 2;
const HighsInt field_2_start = 4;
const HighsInt field_2_width = 8;
const HighsInt field_3_start = 14;
const HighsInt field_3_width = 8;
const HighsInt field_4_start = 24;
const HighsInt field_4_width = 12;
const HighsInt field_5_start = 39;
const HighsInt field_5_width = 8;
const HighsInt field_6_start = 49;
const HighsInt field_6_width = 12;

FilereaderRetcode readMps(
    const HighsLogOptions& log_options, const std::string& filename,
    HighsInt mxNumRow, HighsInt mxNumCol, HighsInt& numRow, HighsInt& numCol,
    ObjSense& objSense, double& objOffset, vector<HighsInt>& Astart,
    vector<HighsInt>& Aindex, vector<double>& Avalue, vector<double>& colCost,
    vector<double>& colLower, vector<double>& colUpper,
    vector<double>& rowLower, vector<double>& rowUpper,
    vector<HighsVarType>& integerColumn, std::string& objective_name,
    vector<std::string>& col_names, vector<std::string>& row_names,
    HighsInt& Qdim, vector<HighsInt>& Qstart, vector<HighsInt>& Qindex,
    vector<double>& Qvalue, HighsInt& cost_row_location, bool& warning_issued,
    const HighsInt keep_n_rows = 0);

HighsStatus writeMps(
    const HighsLogOptions& log_options, const std::string& filename,
    const std::string& model_name, const HighsInt& num_row,
    const HighsInt& num_col, const HighsInt& q_dim, const ObjSense& sense,
    const double& offset, const vector<double>& col_cost,
    const vector<double>& col_lower, const vector<double>& col_upper,
    const vector<double>& row_lower, const vector<double>& row_upper,
    const vector<HighsInt>& a_start, const vector<HighsInt>& a_index,
    const vector<double>& a_value, const vector<HighsInt>& q_start,
    const vector<HighsInt>& q_index, const vector<double>& q_value,
    const vector<HighsVarType>& integrality, const std::string& objective_name,
    const vector<std::string>& col_names, const vector<std::string>& row_names,
    const bool use_free_format = true);

bool load_mpsLine(std::istream& file, HighsVarType& integerVar, HighsInt lmax,
                  char* line, char* flag, double* data);

HighsStatus writeModelAsMps(const HighsOptions& options,
                            const std::string& filename,
                            const HighsModel& model, const bool free = true);

#endif /* IO_HMPSIO_H_ */
