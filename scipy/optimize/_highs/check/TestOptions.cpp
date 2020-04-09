#include <cstdio>

//#include "FilereaderEms.h"
#include "HMPSIO.h"
//#include "HMpsFF.h"
#include "Highs.h"
//#include "HighsIO.h"
//#include "HighsLp.h"
#include "LoadOptions.h"
#include "catch.hpp"

TEST_CASE("internal-options", "[highs_options]") {
  HighsOptions options;
  OptionStatus return_status = checkOptions(options.logfile, options.records);
  REQUIRE(return_status == OptionStatus::OK);

  options.options_file = std::string(HIGHS_DIR) + "/check/sample_options_file";

  bool success = loadOptionsFromFile(options);
  REQUIRE(success == true);
  REQUIRE(options.presolve == on_string);
  REQUIRE(options.small_matrix_value == 0.001);
  REQUIRE(options.mps_parser_type_free);

  reportOptions(stdout, options.records, true);

  return_status = checkOptions(options.logfile, options.records);
  REQUIRE(return_status == OptionStatus::OK);

  // Check setting boolean options
  std::string setting_string = "fixed";
  return_status = setOptionValue(options.logfile, "mps_parser_type_free",
                                 options.records, setting_string);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue(options.logfile, "mps_parser_type_free",
                                 options.records, "fixed");
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue(options.logfile, "mps_parser_type_free",
                                 options.records, "False");
  REQUIRE(return_status == OptionStatus::OK);

  return_status = setOptionValue(options.logfile, "mps_parser_type_free",
                                 options.records, "F");
  REQUIRE(return_status == OptionStatus::OK);

  bool mps_parser_type_free = false;
  return_status = setOptionValue(options.logfile, "mps_parser_type_free",
                                 options.records, mps_parser_type_free);
  REQUIRE(return_status == OptionStatus::OK);

  return_status =
      setOptionValue(options.logfile, "mps_parser_type", options.records, true);
  REQUIRE(return_status == OptionStatus::UNKNOWN_OPTION);

  // Check setting int options

  return_status =
      setOptionValue(options.logfile, "allowed_simplex_matrix_scale_factor",
                     options.records, -1);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status =
      setOptionValue(options.logfile, "allowed_simplex_matrix_scale_factor",
                     options.records, 25);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  std::string allowed_simplex_matrix_scale_factor_string = "1e-7";
  return_status = setOptionValue(
      options.logfile, "allowed_simplex_matrix_scale_factor", options.records,
      allowed_simplex_matrix_scale_factor_string);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status =
      setOptionValue(options.logfile, "allowed_simplex_matrix_scale_factor",
                     options.records, "3.14159");
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  printf("\nAfter setting allowed_simplex_matrix_scale_factor to 1\n");
  reportOptions(stdout, options.records);

  double allowed_simplex_matrix_scale_factor_double = 1e-7;
  return_status = setOptionValue(
      options.logfile, "allowed_simplex_matrix_scale_factor", options.records,
      allowed_simplex_matrix_scale_factor_double);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  int allowed_simplex_matrix_scale_factor = 12;
  return_status =
      setOptionValue(options.logfile, "allowed_simplex_matrix_scale_factor",
                     options.records, allowed_simplex_matrix_scale_factor);
  REQUIRE(return_status == OptionStatus::OK);

  printf("\nAfter testing int options\n");
  reportOptions(stdout, options.records);

  // Check setting double options

  return_status = setOptionValue(options.logfile, "large_matrix_value",
                                 options.records, -1);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue(options.logfile, "large_matrix_value",
                                 options.records, "1");
  REQUIRE(return_status == OptionStatus::OK);

  return_status = setOptionValue(options.logfile, "small_matrix_value",
                                 options.records, -1);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue(options.logfile, "small_matrix_value",
                                 options.records, "1e-6");
  REQUIRE(return_status == OptionStatus::OK);

  double small_matrix_value = 1e-7;
  return_status = setOptionValue(options.logfile, "small_matrix_value",
                                 options.records, small_matrix_value);
  REQUIRE(return_status == OptionStatus::OK);

  // Check setting string options

  return_status = setOptionValue(options.logfile, presolve_string,
                                 options.records, "ml.mps");
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  std::string model_file = "ml.mps";
  return_status = setOptionValue(options.logfile, presolve_string,
                                 options.records, model_file);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status =
      setOptionValue(options.logfile, presolve_string, options.records, "off");
  REQUIRE(return_status == OptionStatus::OK);

  std::string presolve = "choose";
  return_status = setOptionValue(options.logfile, presolve_string,
                                 options.records, presolve);
  REQUIRE(return_status == OptionStatus::OK);

  return_status = setOptionValue(options.logfile, model_file_string,
                                 options.records, model_file);
  REQUIRE(return_status == OptionStatus::OK);

  reportOptions(stdout, options.records);

  bool get_mps_parser_type_free;
  return_status = getOptionValue(options.logfile, "mps_parser_type_free",
                                 options.records, get_mps_parser_type_free);
  REQUIRE(return_status == OptionStatus::OK);
  REQUIRE(get_mps_parser_type_free == false);

  int get_allowed_simplex_matrix_scale_factor;
  return_status =
      getOptionValue(options.logfile, "allowed_simplex_matrix_scale_factor",
                     options.records, get_allowed_simplex_matrix_scale_factor);
  REQUIRE(return_status == OptionStatus::OK);
  REQUIRE(get_allowed_simplex_matrix_scale_factor ==
          allowed_simplex_matrix_scale_factor);

  double get_small_matrix_value;
  return_status = getOptionValue(options.logfile, "small_matrix_value",
                                 options.records, get_small_matrix_value);
  REQUIRE(return_status == OptionStatus::OK);
  REQUIRE(get_small_matrix_value == small_matrix_value);

  std::string get_model_file;
  return_status = getOptionValue(options.logfile, "model_file", options.records,
                                 get_model_file);
  REQUIRE(return_status == OptionStatus::OK);
  REQUIRE(get_model_file == model_file);

  return_status = checkOptions(options.logfile, options.records);
  REQUIRE(return_status == OptionStatus::OK);
}

TEST_CASE("highs-options", "[highs_options]") {
  Highs highs;
  HighsStatus return_status = highs.writeHighsOptions("Highs.set");
  REQUIRE(return_status == HighsStatus::OK);

  // Check setting boolean options
  std::string setting_string = "fixed";
  return_status =
      highs.setHighsOptionValue("mps_parser_type_free", setting_string);
  REQUIRE(return_status == HighsStatus::Error);

  return_status = highs.setHighsOptionValue("mps_parser_type_free", "fixed");
  REQUIRE(return_status == HighsStatus::Error);

  return_status = highs.setHighsOptionValue("mps_parser_type_free", "False");
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.setHighsOptionValue("mps_parser_type_free", "F");
  REQUIRE(return_status == HighsStatus::OK);

  bool mps_parser_type_free = false;
  return_status =
      highs.setHighsOptionValue("mps_parser_type_free", mps_parser_type_free);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.setHighsOptionValue("mps_parser_type", true);
  REQUIRE(return_status == HighsStatus::Error);

  // Check setting int options

  return_status =
      highs.setHighsOptionValue("allowed_simplex_matrix_scale_factor", -1);
  REQUIRE(return_status == HighsStatus::Error);

  return_status =
      highs.setHighsOptionValue("allowed_simplex_matrix_scale_factor", 25);
  REQUIRE(return_status == HighsStatus::Error);

  std::string allowed_simplex_matrix_scale_factor_string = "1e-7";
  return_status =
      highs.setHighsOptionValue("allowed_simplex_matrix_scale_factor",
                                allowed_simplex_matrix_scale_factor_string);
  REQUIRE(return_status == HighsStatus::Error);

  return_status = highs.setHighsOptionValue(
      "allowed_simplex_matrix_scale_factor", "3.14159");
  REQUIRE(return_status == HighsStatus::Error);

  printf("\nAfter setting allowed_simplex_matrix_scale_factor to 1\n");
  return_status = highs.writeHighsOptions("Highs.set");
  REQUIRE(return_status == HighsStatus::OK);

  double allowed_simplex_matrix_scale_factor_double = 1e-7;
  return_status =
      highs.setHighsOptionValue("allowed_simplex_matrix_scale_factor",
                                allowed_simplex_matrix_scale_factor_double);
  REQUIRE(return_status == HighsStatus::Error);

  int allowed_simplex_matrix_scale_factor = 12;
  return_status =
      highs.setHighsOptionValue("allowed_simplex_matrix_scale_factor",
                                allowed_simplex_matrix_scale_factor);
  REQUIRE(return_status == HighsStatus::OK);

  printf("\nAfter testing int options\n");
  return_status = highs.writeHighsOptions("Highs.set");
  REQUIRE(return_status == HighsStatus::OK);

  // Check setting double options

  return_status = highs.setHighsOptionValue("large_matrix_value", -1);
  REQUIRE(return_status == HighsStatus::Error);

  return_status = highs.setHighsOptionValue("large_matrix_value", "1");
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.setHighsOptionValue("small_matrix_value", -1);
  REQUIRE(return_status == HighsStatus::Error);

  return_status = highs.setHighsOptionValue("small_matrix_value", "1e-6");
  REQUIRE(return_status == HighsStatus::OK);

  double small_matrix_value = 1e-7;
  return_status =
      highs.setHighsOptionValue("small_matrix_value", small_matrix_value);
  REQUIRE(return_status == HighsStatus::OK);

  // Check setting string options

  return_status = highs.setHighsOptionValue(presolve_string, "ml.mps");
  REQUIRE(return_status == HighsStatus::Error);

  std::string model_file = "ml.mps";
  return_status = highs.setHighsOptionValue(presolve_string, model_file);
  REQUIRE(return_status == HighsStatus::Error);

  return_status = highs.setHighsOptionValue(presolve_string, "off");
  REQUIRE(return_status == HighsStatus::OK);

  std::string presolve = "choose";
  return_status = highs.setHighsOptionValue(presolve_string, presolve);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.setHighsOptionValue(model_file_string, model_file);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.writeHighsOptions("Highs.set");
  REQUIRE(return_status == HighsStatus::OK);

  bool get_mps_parser_type_free;
  return_status = highs.getHighsOptionValue("mps_parser_type_free",
                                            get_mps_parser_type_free);
  REQUIRE(return_status == HighsStatus::OK);
  REQUIRE(get_mps_parser_type_free == false);

  int get_allowed_simplex_matrix_scale_factor;
  return_status =
      highs.getHighsOptionValue("allowed_simplex_matrix_scale_factor",
                                get_allowed_simplex_matrix_scale_factor);
  REQUIRE(return_status == HighsStatus::OK);
  REQUIRE(get_allowed_simplex_matrix_scale_factor ==
          allowed_simplex_matrix_scale_factor);

  double get_small_matrix_value;
  return_status =
      highs.getHighsOptionValue("small_matrix_value", get_small_matrix_value);
  REQUIRE(return_status == HighsStatus::OK);
  REQUIRE(get_small_matrix_value == small_matrix_value);

  std::string get_model_file;
  return_status = highs.getHighsOptionValue("model_file", get_model_file);
  REQUIRE(return_status == HighsStatus::OK);
  REQUIRE(get_model_file == model_file);

  HighsOptions options = highs.getHighsOptions();
  REQUIRE(options.model_file == model_file);
  REQUIRE(options.small_matrix_value == small_matrix_value);
  REQUIRE(options.allowed_simplex_matrix_scale_factor ==
          allowed_simplex_matrix_scale_factor);
  REQUIRE(options.mps_parser_type_free == mps_parser_type_free);
}
