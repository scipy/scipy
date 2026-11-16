/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsDeprecated.cpp
 * @brief
 */
#include "HConfig.h"
#include "Highs.h"

HighsStatus Highs::setLogCallback(void (*user_log_callback)(HighsLogType,
                                                            const char*, void*),
                                  void* user_log_callback_data) {
  deprecationMessage("setLogCallback", "setCallback");
  options_.log_options.user_log_callback = user_log_callback;
  options_.log_options.user_log_callback_data = user_log_callback_data;
  return HighsStatus::kOk;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const bool value) {
  deprecationMessage("setHighsOptionValue", "setOptionValue");
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const HighsInt value) {
  deprecationMessage("setHighsOptionValue", "setOptionValue");
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const double value) {
  deprecationMessage("setHighsOptionValue", "setOptionValue");
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const std::string& value) {
  deprecationMessage("setHighsOptionValue", "setOptionValue");
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const char* value) {
  deprecationMessage("setHighsOptionValue", "setOptionValue");
  return setOptionValue(option, value);
}

HighsStatus Highs::readHighsOptions(const std::string& filename) {
  deprecationMessage("readHighsOptions", "readOptions");
  return readOptions(filename);
}

HighsStatus Highs::passHighsOptions(const HighsOptions& options) {
  deprecationMessage("passHighsOptions", "passOptions");
  return passOptions(options);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option, bool& value) {
  deprecationMessage("getHighsOptionValue", "getOptionValue");
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       HighsInt& value) {
  deprecationMessage("getHighsOptionValue", "getOptionValue");
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       double& value) {
  deprecationMessage("getHighsOptionValue", "getOptionValue");
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       std::string& value) {
  deprecationMessage("getHighsOptionValue", "getOptionValue");
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionType(const std::string& option,
                                      HighsOptionType& type) {
  deprecationMessage("getHighsOptionType", "getOptionType");
  return getOptionType(option, type);
}

HighsStatus Highs::resetHighsOptions() {
  deprecationMessage("resetHighsOptions", "resetOptions");
  return resetOptions();
}

HighsStatus Highs::writeHighsOptions(
    const std::string& filename, const bool report_only_non_default_values) {
  deprecationMessage("writeHighsOptions", "writeOptions");
  return writeOptions(filename, report_only_non_default_values);
}

const HighsOptions& Highs::getHighsOptions() const {
  deprecationMessage("getHighsOptions", "getOptions");
  return getOptions();
}

HighsStatus Highs::setHighsLogfile(FILE* logfile) {
  deprecationMessage("setHighsLogfile", "None");
  options_.output_flag = false;
  return HighsStatus::kOk;
}

HighsStatus Highs::setHighsOutput(FILE* output) {
  deprecationMessage("setHighsOutput", "None");
  options_.output_flag = false;
  return HighsStatus::kOk;
}

const HighsInfo& Highs::getHighsInfo() const {
  deprecationMessage("getHighsInfo", "getInfo");
  return getInfo();
}

HighsStatus Highs::getHighsInfoValue(const std::string& info, HighsInt& value) {
  deprecationMessage("getHighsInfoValue", "getInfoValue");
  return getInfoValue(info, value);
}

HighsStatus Highs::getHighsInfoValue(const std::string& info,
                                     double& value) const {
  deprecationMessage("getHighsInfoValue", "getInfoValue");
  return getInfoValue(info, value);
}

HighsStatus Highs::writeHighsInfo(const std::string& filename) {
  deprecationMessage("writeHighsInfo", "writeInfo");
  return writeInfo(filename);
}

double Highs::getHighsInfinity() {
  deprecationMessage("getHighsInfinity", "getInfinity");
  return getInfinity();
}

double Highs::getHighsRunTime() {
  deprecationMessage("getHighsRunTime", "getRunTime");
  return getRunTime();
}

#if 0
HighsStatus Highs::writeSolution(const std::string& filename,
                                const bool pretty) const {
  deprecationMessage("writeSolution(filename, pretty)", "writeSolution(filename, style)");
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  FILE* file;
  HighsFileType file_type;
  call_status = openWriteFile(filename, "writeSolution", file, file_type);
  return_status =
      interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::kError) return return_status;
  HighsInt style;
  if (pretty) {
    style = kSolutionStylePretty;
  } else {
    style = kSolutionStyleRaw;
  }
  writeSolutionFile(file, options_,
		    model_, basis_, solution_, info_, model_status_,
                    style);
  if (file != stdout) fclose(file);
  return HighsStatus::kOk;
}
#endif

const HighsModelStatus& Highs::getModelStatus(const bool) const {
  deprecationMessage("getModelStatus(const bool scaled_model)",
                     "getModelStatus()");
  return model_status_;
}
