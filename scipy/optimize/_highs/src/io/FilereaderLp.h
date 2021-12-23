/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#ifndef IO_FILEREADER_LP_H_
#define IO_FILEREADER_LP_H_

#include <list>

#include "io/Filereader.h"
#include "io/HighsIO.h"

#define BUFFERSIZE 561
#define LP_MAX_LINE_LENGTH 560
#define LP_MAX_NAME_LENGTH 255

#define LP_COMMENT_FILESTART ("File written by Highs .lp filereader")

class FilereaderLp : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const HighsOptions& options,
                                      HighsLp& model);

  HighsStatus writeModelToFile(const HighsOptions& options,
                               const std::string filename, HighsLp& model);

 private:
  // functions to write files
  int linelength;
  void writeToFile(FILE* file, const char* format, ...);
  void writeToFileLineend(FILE* file);
};

#endif
