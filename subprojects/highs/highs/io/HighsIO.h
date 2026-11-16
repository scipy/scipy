/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HighsIO.h
 * @brief IO methods for HiGHS - currently just print/log messages
 */
#ifndef HIGHS_IO_H
#define HIGHS_IO_H

#include <array>
#include <iostream>

#include "lp_data/HighsCallback.h"

class HighsOptions;

const HighsInt kIoBufferSize = 1024;  // 65536;

enum class HighsFileType { kMinimal = 0, kFull, kMps, kLp, kMd };

/**
 * @brief IO methods for HiGHS - currently just print/log messages
 */
const char* const HighsLogTypeTag[] = {"", "",          "",
                                       "", "WARNING: ", "ERROR:   "};
enum LogDevLevel {
  kHighsLogDevLevelMin = 0,
  kHighsLogDevLevelNone = kHighsLogDevLevelMin,  // 0
  kHighsLogDevLevelInfo,                         // 1
  kHighsLogDevLevelDetailed,                     // 2
  kHighsLogDevLevelVerbose,                      // 3
  kHighsLogDevLevelMax = kHighsLogDevLevelVerbose
};

struct HighsLogOptions {
  FILE* log_stream;
  bool* output_flag;
  bool* log_to_console;
  HighsInt* log_dev_level;
  void (*user_log_callback)(HighsLogType, const char*, void*);
  void* user_log_callback_data;
  HighsCallbackFunctionType user_callback;
  void* user_callback_data;
  bool user_callback_active;
  void clear();
  HighsLogOptions()
      : log_stream(nullptr),
        output_flag(nullptr),
        log_to_console(nullptr),
        log_dev_level(nullptr),
        user_log_callback(nullptr),
        user_log_callback_data(nullptr),
        user_callback_data(nullptr),
        user_callback_active(false){};
};

/**
 * @brief Write the HiGHS version and copyright statement
 */
void highsLogHeader(const HighsLogOptions& log_options, const bool log_githash);

/**
 * @brief Convert a double number to a string using given tolerance
 */
std::array<char, 32> highsDoubleToString(const double val,
                                         const double tolerance);

/**
 * @brief Convert a MIP relative gap number to a string
 */
std::array<char, 32> highsGapToString(const double gap);

/**
 * @brief For _single-line_ user logging with message type notification
 */
// Printing format: must contain exactly one "\n" at end of format
void highsLogUser(const HighsLogOptions& log_options_, const HighsLogType type,
                  const char* format, ...);

/**
 * @brief For development logging
 */
void highsLogDev(const HighsLogOptions& log_options_, const HighsLogType type,
                 const char* format, ...);

/**
 * @brief Gets filename extension
 */
std::string getFilenameExt(const std::string& filename);

/**
 * @brief Replaces fprintf(file,... so that when file=stdout highsLogUser is
 * used
 */
// Printing format: must contain exactly one "\n" at end of format
void highsFprintfString(FILE* file, const HighsLogOptions& log_options_,
                        const std::string& s);

/**
 * @brief For development logging when true log_options may not be available -
 * indicated by null pointer
 */
void highsReportDevInfo(const HighsLogOptions* log_options,
                        const std::string& line);

void highsOpenLogFile(HighsOptions& options, const std::string& log_file);

void highsReportLogOptions(const HighsLogOptions& log_options_);
void highsPrintStringAsHeading(const std::string& text);

std::string highsFormatToString(const char* format, ...);

const std::string highsBoolToString(const bool b,
                                    const HighsInt field_width = 2);

const std::string highsIntToPlural(const HighsInt i, const bool y = false);
const std::string highsIntToBe(const HighsInt i);
const std::string highsInsertMdEscapes(const std::string& from_string);
const std::string highsInsertMdId(const std::string& from_string);
const std::string highsTimeToString(const double time);
const std::string highsTimeSecondToString(const double time);
#endif
