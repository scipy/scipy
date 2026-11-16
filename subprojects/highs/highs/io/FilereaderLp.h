/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief
 */

#ifndef IO_FILEREADER_LP_H_
#define IO_FILEREADER_LP_H_

#include <list>

#include "io/Filereader.h"
#include "io/HighsIO.h"

#define BUFFERSIZE 561
#define LP_MAX_LINE_LENGTH 560
#define LP_MAX_NAME_LENGTH 255

#define LP_COMMENT_FILESTART ("File written by HiGHS .lp file handler")

class FilereaderLp : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const HighsOptions& options,
                                      const std::string filename,
                                      HighsModel& model);

  HighsStatus writeModelToFile(const HighsOptions& options,
                               const std::string filename,
                               const HighsModel& model);

 private:
  // functions to write files
  HighsInt linelength;
  void writeToFile(FILE* file, const char* format, ...);
  void writeToFileLineEnd(FILE* file);
  void writeToFileValue(FILE* file, const double value,
                        const bool force_plus = true);
  void writeToFileVar(FILE* file, const std::string& var_name);
  void writeToFileMatrixRow(FILE* file, const HighsInt iRow,
                            const HighsSparseMatrix& ar_matrix,
                            const std::vector<string>& col_names);
};

#endif
