/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/LoadOptions.h
 * @brief
 */

#ifndef IO_LOAD_OPTIONS_H_
#define IO_LOAD_OPTIONS_H_

#include "lp_data/HighsOptions.h"

enum class HighsLoadOptionsStatus { kError = -1, kOk = 0, kEmpty = 1 };

// For extended options to be parsed from filename
HighsLoadOptionsStatus loadOptionsFromFile(
    const HighsLogOptions& report_log_options, HighsOptions& options,
    const std::string& filename);

#endif
