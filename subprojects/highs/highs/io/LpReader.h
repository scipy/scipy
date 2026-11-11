/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/LpReader.h
 * @brief Reader for the CPLEX LP file format
 */
#ifndef IO_LP_READER_H_
#define IO_LP_READER_H_

#include <string>

#include "io/Filereader.h"
#include "io/HighsIO.h"
#include "model/HighsModel.h"

// Read the LP file `filename`, which may be compressed, into `model`, which
// must be empty. The file is read in chunks, and the model is built as the
// file is parsed.
//
// Syntax errors are reported through `log_options` as diagnostics that give
// the location of the error in the file. Returns:
//
//  - FilereaderRetcode::kOk if the file was read without problems
//  - FilereaderRetcode::kWarning if the file was read, but something in it was
//    ignored or modified
//  - FilereaderRetcode::kFileNotFound if the file cannot be opened
//  - FilereaderRetcode::kParserError if the file is not a valid LP file, or it
//    contains something HiGHS cannot handle, in which case `model` is
//    incomplete
FilereaderRetcode readLpFile(const HighsLogOptions& log_options,
                             const std::string& filename, HighsModel& model);

#endif
